import os, sys, subprocess, multiprocessing, concurrent.futures, time, random, copy, itertools
import functools
from tqdm.auto import tqdm

"""
Polyalign
Polyalign is a a wrapper for bwa mem (and bwa-mem2) with improved read pairing.
Left and right reads are aligned separately, and all pairings with correct insert size are returned.
"""

class BwaMem:
    """
    Object to handle bwa index indexing and bwa mem alignment.
    Handles either input fastq files or lists of lines, outputing sam output lines (as split lists) as an iterator.
    """
    def __init__(self, subject_path, options={}, program="bwa", do_indexing=True, cpus=None):
        self.subject_path = subject_path
        self.options = options
        programs = ["bwa", "bwa-mem2"]
        self.program = program
        if program not in programs:
            raise ValueError("Program must be one of: "+", ".join(programs))
        self.cpus = cpus
        if self.cpus is None:
            self.cpus = multiprocessing.cpu_count() // 2
        self.cpus = str(self.cpus)
        self.index_suffices = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        if do_indexing:
            self.makeIndex()
    
    def checkBwa(self):
        """
        Calls the program to check that it exists.
        """
        try:
            subprocess.run([self.program, "mem", "help"], stdout=subprocess.PIPE)
        except OSError:
            print("[PA::BwaMem]", "Error: " + self.program + " not found")

    def makeIndex(self):
        """
        Index subject_path fasta file using bwa index.
        """
        index_exists = True
        for suffix in self.index_suffices:
            if not os.path.exists(self.subject_path + suffix):
                index_exists = False
                print("[PA::BwaMem]", "Index files not found,", self.subject_path + suffix, "missing")
                break
        if not index_exists:
            subprocess.run([self.program, "index", self.subject_path])
    
    def removeIndex(self):
        """
        Removes index files for subject_path.
        """
        for suffix in self.index_suffices:
            try:
                os.remove(self.subject_path + suffix)
            except FileNotFoundError:
                pass
    
    class SamLineIterator:
        """
        Iterator object for returning sam lines as lists, with simple caching of current value.
        """
        def __init__(self, iterator):
            self.iterator = iterator
            self.current = None
            self.current_id = None
            self.previous = None
            self.previous_id = None
        
        def __iter__(self):
            return self
        
        def __next__(self):
            try:
                self.previous = self.current
                self.previous_id = self.current_id
                self.current = self.convert(next(self.iterator))
                if self.current == [""]:
                    raise StopIteration
                if self.current[0][0] != "@":
                    self.current_id = self.current[0]
                else:
                    self.current_id = None
            except StopIteration:
                self.current = None
                self.current_id = None
            finally:
                return self.current
        
        def convert(self, entry):
            """
            Convert bytestream to line split by tab delimitation.
            """
            if entry is None:
                return None
            return entry.decode("utf-8").strip().split("\t")
        
        def readLine(self):
            """
            Read and return the next line
            """
            return next(self)
        
        def readHeader(self):
            """
            Read and return all header lines as a list of lists
            Does not read first line, should be run after next(<SamLineIterator instance>) or readLine()
            """
            line = self.current
            if line is None:
                return None
            lines = [line.copy()]
            while next(self) is not None and self.current[0][0] == "@":
                lines.append(self.current.copy())
            return lines
        
        def nextRead(self):
            """
            Return list of all sam lines for next read id as a list of lists
            Does not check for header lines, should be run after readHeader()
            """
            line = self.current
            if line is None:
                return None
            lines = [line.copy()]
            while next(self) is not None and lines[0][0]==self.current_id:
                lines.append(self.current.copy())
            return lines
    
    def mem(self, reads_1_path, reads_2_path=None):
        """
        Do alignment using bwa mem of reads_1_path (and optionally reads_2_path).
        """
        command = [self.program, "mem", "-t", self.cpus, self.subject_path]
        for option in self.options:
            if self.options[option] is not None:
                command += ["-" + option, self.options[option]]
            else:
                command += ["-" + option]
        command += [reads_1_path]
        if reads_2_path is not None:
            command += [reads_2_path]
        print("[PA::BwaMem]", "bwa mem command:", " ".join(command))
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        return self.SamLineIterator(iter(process.stdout.readline, ""))
    
    def memString(self, reads_1_lines, reads_2_lines=None, temp_directory_path=None, options={}):
        """
        Do alignment using bwa mem of list of data lines reads_1_lines (and optionally reads_2_lines)
        Construct temporary input files to do searching in temp_path.
        """
        if temp_directory_path is None:
            temp_directory_path = os.path.dirname(self.subject_path)
        files = {"1": {"data": reads_1_lines}}
        if reads_2_lines is not None:
            files = {"2": {"data": reads_2_lines}}
        for key in files:
            files[key]["path"] = os.path.join(temp_directory_path, "reads_temp_"+key+".fastq")
            with open(files[key]["path"], "w") as file:
                for line in files[key]["data"]:
                    file.write(line + "\n")
        if reads_2_lines is not None:
            self.mem(files["1"]["path"], files["2"]["path"], options)
        else:
            self.mem(files["1"]["path"], options)

class Header:
    """
    Object to contain sam header data
    """
    def __init__(self, sam_line):
        self.line = sam_line

    def __str__(self):
        """
        Return sam line as a string
        """
        return "\t".join(self.columns())

    def columns(self):
        """
        Sam line as a list of strings for each column
        """
        return self.line

class Alignment:
    """
    Object to contain read data as derived from sam-formatted line
    https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    def __init__(self, sam_line):
        self.sam_line = sam_line
        self.flag_indices = ["primary", "properly_aligned", "unmapped", "mate_unmapped", "reverse", "mate_reverse", "first", "second", "secondary", "qc_fail", "duplicate", "supplementary"]
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """
        Return sam line as a string
        """
        return "\t".join(self.columns())
    
    @functools.cached_property
    def qname(self):
        return self.sam_line[0]
    
    @functools.cached_property
    def flag(self):
        return int(self.sam_line[1])
    
    @functools.cached_property
    def rname(self):
        return self.sam_line[2]
    
    @functools.cached_property
    def pos(self):
        return int(self.sam_line[3])
    
    @functools.cached_property
    def mapq(self):
        return int(self.sam_line[4])
    
    @functools.cached_property
    def cigar(self):
        return self.sam_line[5]
    
    @functools.cached_property
    def rnext(self):
        return self.sam_line[6]
    
    @functools.cached_property
    def pnext(self):
        return int(self.sam_line[7])
    
    @functools.cached_property
    def tlen(self):
        return int(self.sam_line[8])
    
    @functools.cached_property
    def seq(self):
        return self.sam_line[9]
    
    @functools.cached_property
    def qual(self):
        return self.sam_line[10]
    
    @functools.cached_property
    def optional(self):
        return self.sam_line[11:]

    def columns(self):
        """
        Sam line as a list of strings for each column
        """
        return [self.qname, str(self.flag), self.rname, str(self.pos), str(self.mapq), self.cigar, self.rnext, str(self.pnext), str(self.tlen), self.seq, self.qual] + self.optional

    def position(self):
        """
        0 indexed start, end position tuple
        start < end for forward, start > end for reverse
        """
        #start, end = self.pos - 1, self.pos - 1 + self.tlen # tlen is inconsistently set(?)
        start, end = self.pos - 1, self.pos - 1 + self.cigarTemplateLength()
        if self.orientation() == "f":
            return start, end
        return end, start
    
    def getFlags(self):
        """
        Flags as a dict of boolean values
        """
        flags = {}
        for index, flag in enumerate(self.flag_indices):
            flags[flag] = self.flag & 2**index != 0
        return flags

    def setFlags(self, flags):
        """
        Set flags from a dict of boolean values
        """
        for index, flag in enumerate(self.flag_indices):
            if flags[flag]:
                self.flag |= 2**index
            else:
                self.flag &= ~2**index
    
    def setFlag(self, flag, value):
        """
        Set a single flag by name
        """
        index = self.flag_indices.index(flag)
        if value:
            self.flag |= 2**index
        else:
            self.flag &= ~2**index
    
    def setOptional(self, tag, value, typechar="Z"):
        """
        Add an optional field to the end of the line
        """
        for entry in self.optional:
            if entry.split(":")[0] == tag:
                self.optional.remove(entry)
        self.optional.append(tag + ":" + typechar + ":" + value)

    def leftmostPosition(self):
        return min(self.position())
        
    def rightmostPosition(self):
        return max(self.position())
    
    def orientation(self):
        """
        f for forward, r for reverse
        """
        if self.flag & 0x10 == 0:
            return "f"
        return "r"
    
    def aligned(self):
        if self.rname == "*":
            return False
        return True
    
    def cigarOperations(self):
        """
        List of cigar operations as tuples of (operation, length)
        """
        cigar = self.cigar
        operations = []
        current_length = ""
        for char in cigar:
            if char.isdigit():
                current_length += char
            else:
                operations.append((char, int(current_length)))
                current_length = ""
        return operations

    def cigarTemplateLength(self):
        """
        Length of the template from the cigar operations
        """
        return sum([length for operation, length in self.cigarOperations() if operation in "MDN"])

    def cigarReadLength(self):
        """
        Length of the query from the cigar operations
        """
        return sum([length for operation, length in self.cigarOperations() if operation in "MIS"])

class AlignmentPair:
    """
    Properties about a mate pair alignment from two sam-formatted lines
    """
    def __init__(self, read_1, read_2, correct_orientation=None, insert_size_range=None):
        self.read_1 = read_1
        self.read_2 = read_2
        self.correct_orientation = correct_orientation
        self.insert_size_range = insert_size_range
    
    def sameReference(self):
        """
        Aligned to the same reference
        """
        if self.read_1.rname == "*" or self.read_2.rname == "*":
            # one or both not aligned
            return False
        if self.read_1.rname == self.read_2.rname:
            # matching reference name
            return True
        return False
    
    def insertSize(self):
        """
        Insert size, from alignment coordinates
        """
        if not self.sameReference():
            return None
        return max(self.read_1.rightmostPosition(), self.read_2.rightmostPosition()) - min(self.read_1.leftmostPosition(), self.read_2.leftmostPosition())
    
    def readLength(self):
        """
        List of both read lengths
        """
        return [len(self.read_1.seq), len(self.read_2.seq)]

    def orientation(self):
        """
        Pair orientation: fr, ff, rr or rf
        """
        if self.read_1.rname == "*" or self.read_2.rname == "*":
            # one or both not aligned
            return None
        if self.read_1.leftmostPosition() < self.read_2.leftmostPosition():
            return self.read_1.orientation() + self.read_2.orientation()
        return self.read_2.orientation() + self.read_1.orientation()
    
    def correctInsertSize(self):
        """
        Is the insert size within the expected range
        """
        if self.insert_size_range is None:
            return None
        if self.insertSize() is None:
            return False
        if self.insert_size_range[0] < self.insertSize() < self.insert_size_range[1]:
            return True
        return False
    
    def correctOrientation(self):
        """
        Is the mate pair oriented correctly
        """
        if self.correct_orientation is None:
            return None
        if self.orientation() == self.correct_orientation:
            return True
        return False

class AlignmentPairFilter:
    def __init__(self, read_length=None, read_orientation=None, insert_size_range=None, retain_unmapped=None, output_type=None):
        self.read_length = read_length
        self.read_orientation = read_orientation
        self.insert_size_range = insert_size_range
        self.retain_unmapped = retain_unmapped
        self.output_type = output_type
    
    def isSinglePair(self, read_1, read_2):
        """
        Check if read_1 and read_2 have length 1, ie. single alignment of paired reads
        """
        if len(read_1) == 1 and len(read_2) == 1:
            return True
        return False
    
    def isUnalignedPair(self, read_1, read_2):
        """
        Check if read_1 and read_2 are both unaligned, this can only occur when both are single alignments
        """
        if self.isSinglePair(read_1, read_2) and not read_1[0].aligned() and not read_2[0].aligned():
            return True
        return False
    
    def isOneUnalignedPair(self, read_1, read_2):
        """
        Check if one of read_1 or read_2 is unaligned, the unaligned read must be the first entry in either read_1 or read_2
        """
        if not read_1[0].aligned() or not read_2[0].aligned():
            return True
        return False
    
    def isGoodSinglePair(self, read_1, read_2):
        """
        Check if read_1 and read_2 are a good single pair, ie. align to the same reference, with the correct insert size, in the correct orientation
        """
        if self.isSinglePair(read_1, read_2):
            read_pair = AlignmentPair(read_1[0], read_2[0], self.read_orientation, self.insert_size_range)
            if read_pair.sameReference() and read_pair.correctInsertSize() and read_pair.correctOrientation():
                return True
        return False
    
    def findGoodPairsNaive(self, read_1, read_2):
        """
        Find all good pairs of read_1 and read_2 alignments, pairs which align to the same reference, with the correct insert size, in the correct orientation
        Naive implementation which exhaustively tests all combinations of read_1 and read_2, checking for same reference, correct insert size and correct orientation
        Scales by len(read_1) * len(read_2)
        """
        good_1, good_2 = [False] * len(read_1), [False] * len(read_2)
        for index_1, alignment_1 in enumerate(read_1):
            for index_2, alignment_2 in enumerate(read_2):
                if self.isGoodSinglePair([alignment_1], [alignment_2]):
                    good_1[index_1] = True
                    good_2[index_2] = True
        result_1, result_2 = [read_1[index] for index, value in enumerate(good_1) if value], [read_2[index] for index, value in enumerate(good_2) if value]
        result_1[0].seq, result_1[0].qual = read_1[0].seq, read_1[0].qual
        result_2[0].seq, result_2[0].qual = read_2[0].seq, read_2[0].qual
        return result_1, result_2
    
    def findGoodPairs(self, read_1, read_2):
        """
        Find all good pairs of read_1 and read_2 alignments, pairs which align to the same reference, with the correct insert size, in the correct orientation
        Smartish implementation which should be much faster for repetitive genomes, by splitting reads by reference name and sorting by position
        For each read_1, only checks read_2 pre-split to the correct reference and carefully stepping list indices to check within the potential good insert size range
        """
        def split_reads_by_rname(read):
            dict = {}
            for alignment in read:
                if alignment.rname != "*":
                    if alignment.rname not in dict:
                        dict[alignment.rname] = []
                    dict[alignment.rname].append(alignment)
            return dict
        
        # setup output
        result_1, result_2 = [], []
        # check for no possible pairs
        if len(read_1) == 0 or len(read_2) == 0:
            return result_1, result_2
        # good pairs must be aligned to the same reference, so start by splitting by rname
        split_read_1 = split_reads_by_rname(read_1)
        split_read_2 = split_reads_by_rname(read_2)
        # search per reference
        for rname in split_read_1:
            # good pairs only possible on same reference
            if rname not in split_read_2:
                continue
            # sort reads by increasing position on rname
            split_read_1[rname].sort(key=lambda x: x.pos)
            split_read_2[rname].sort(key=lambda x: x.pos)
            # do smartish iteration through reads to find good pairs
            good_1, good_2 = [False] * len(split_read_1[rname]), [False] * len(split_read_2[rname])
            index_1, index_2 = 0, 0
            while index_1 < len(split_read_1[rname]) and index_2 < len(split_read_2[rname]):
                # first, increase index_2 until that read_2 is could be within the insert size range
                while index_2 < len(split_read_2[rname]) and split_read_2[rname][index_2].pos + self.insert_size_range[1] + self.read_length < split_read_1[rname][index_1].pos:
                    index_2 += 1
                # if index_2 is now out of range, break
                if index_2 >= len(split_read_2[rname]):
                    break
                # for each read_2 within the potential insert size range, check if it would be a good pair
                cur_index_2 = index_2
                while cur_index_2 < len(split_read_2[rname]) and split_read_2[rname][cur_index_2].pos - self.insert_size_range[1] - self.read_length < split_read_1[rname][index_1].pos:
                    if self.isGoodSinglePair([split_read_1[rname][index_1]], [split_read_2[rname][cur_index_2]]):
                        good_1[index_1], good_2[cur_index_2] = True, True
                    cur_index_2 += 1
                # increase index_1
                index_1 += 1
                # increase index_2 to the next read_2 not yet marked as part of a good pair
                while index_2 < len(split_read_2[rname]) and good_2[index_2]:
                    index_2 += 1
            # add good pairs to output
            result_1 += [split_read_1[rname][index] for index, value in enumerate(good_1) if value]
            result_2 += [split_read_2[rname][index] for index, value in enumerate(good_2) if value]
        # ensure that sequence and quality are set for first retained read
        if len(result_1) > 0:
            result_1[0].seq, result_1[0].qual = read_1[0].seq, read_1[0].qual
        if len(result_2) > 0:
            result_2[0].seq, result_2[0].qual = read_2[0].seq, read_2[0].qual
        return result_1, result_2

    def analyseReadPairing(self, next_read_1, next_read_2):
        """
        Read input sams, check for read alignment and return correctly paired reads.
        Aims to copy Polypolish filtering: https://github.com/rrwick/Polypolish/blob/main/src/filter.rs alignment_pass_qc()
        """
        # no next read, return None
        if next_read_1 is None or next_read_2 is None:
            status = None
            return None, None, status
        # TODO: More consistent handling of read retention/discard, eg. just pop from list if not retained, and add support for setting ZP:Z:fail flag
        # check for read name mismatch
        if next_read_1[0].qname != next_read_2[0].qname:
            print("[PA::ARP]", "Error: Read names do not match")
        # check for no alignments
        if self.isUnalignedPair(next_read_1, next_read_2):
            # For pairs where both have zero alignments, retain sam lines, but no way to contribute to polishing.
            status = "both unaligned"
            next_read_1[0].setOptional("ZP", "fail", "Z")
            next_read_2[0].setOptional("ZP", "fail", "Z")
            if self.retain_unmapped:
                return next_read_1, next_read_2, status
            return [], [], status
        if self.isOneUnalignedPair(next_read_1, next_read_2):
            # For pairs with exactly one alignment of each read, keep sam lines, ie. treat unique paired alignments as correct.
            status = "single unaligned"
            return next_read_1, next_read_2, status
        if self.isSinglePair(next_read_1, next_read_2):
            # For pairs where one has zero alignments, keep sam lines, ie. treat single alignments as correct.
            status = "both unique aligned"
            return next_read_1, next_read_2, status
        # For pairs with multiple alignments, keep sam lines if it pairs (same ref seq, correct ori, good insert) with any of its pair's alignments.
        good_1, good_2 = self.findGoodPairs(next_read_1, next_read_2)
        if self.retain_unmapped:
            # Keep all read alignments, but mark as failed
            for read_1 in next_read_1:
                if read_1 not in good_1:
                    read_1.setOptional("ZP", "fail", "Z")
            for read_2 in next_read_2:
                if read_2 not in good_2:
                    read_2.setOptional("ZP", "fail", "Z")
            if len(good_1) == 0 or len(good_2) == 0:
                status = "no good pairs from multiple"
            else:
                status = "good pairs from multiple"
            return next_read_1, next_read_2, status
        else:
            # keep only good read alignments
            if len(good_1) == 0 or len(good_2) == 0:
                status = "no good pairs from multiple"
                return [], [], status
            for index, good in enumerate(good_1):
                # add sequence and quality to first good reads, replace with "*" for the rest
                if index == 0:
                    good.seq = next_read_1[0].seq
                    good.qual = next_read_1[0].qual
                else:
                    good.seq = "*"
                    good.qual = "*"
            for index, good in enumerate(good_2):
                if index == 0:
                    good.seq = next_read_2[0].seq
                    good.qual = next_read_1[0].qual
                else:
                    good.seq = "*"
                    good.qual = "*"
            status = "good pairs from multiple"
            return good_1, good_2, status
        
    def alignmentOutputPreprocess(self, lines_1, lines_2):
        """
        Prepare alignment lines for output based on output type
        """
        if self.output_type == "filtered":
            output = [[], []]
            for l in range(len(lines_1)):
                output[0] += [str(line) for line in lines_1[l]]
                output[1] += [str(line) for line in lines_2[l]]
            return output
        if self.output_type == "paired":
            output = [[]]
            for l in range(len(lines_1)):
                # for each read pair
                if len(lines_1[l]) > 0 and len(lines_2[l]) > 0:
                    # pick a random pair if multiple alignments, ensuring sequence data is included
                    if len(lines_1[l]) > 1 or len(lines_2[l]) > 1:
                        out_1, out_2 = random.choice(lines_1[l]), random.choice(lines_2[l])
                        out_1.seq, out_1.qual = lines_1[l][0].seq, lines_1[l][0].qual
                        out_2.seq, out_2.qual = lines_2[l][0].seq, lines_2[l][0].qual
                        lines_1[l], lines_2[l] = [out_1], [out_2]
                    # set pairing information
                    lines_1[l][0].rnext, lines_2[l][0].rnext = lines_2[l][0].rname, lines_1[l][0].rname
                    lines_1[l][0].pnext, lines_2[l][0].pnext = lines_2[l][0].pos, lines_1[l][0].pos
                    output[0].append(str(lines_1[l][0]))
                    output[0].append(str(lines_2[l][0]))
            return output
        if self.output_type == "splitfiltered" or self.output_type == "filtered":
            output = [{}, {}]
            lines = [lines_1, lines_2]
            for dir in range(len(lines)):
                # for left and right read
                for l in range(len(lines[dir])):
                    # for each read
                    out_rnames = [] # each reference aligned to, so far
                    for a in range(len(lines[dir][l])):
                        # for each alignment
                        rname = lines[dir][l][a].rname
                        if rname != "*": # only write to buffer if reference is not *, ie. unaligned sequcence
                            if rname not in out_rnames:
                                # first alignment to this reference, so infill sequence data
                                out_rnames.append(rname)
                                lines[dir][l][a].seq, lines[dir][l][a].qual = lines[dir][l][0].seq, lines[dir][l][0].qual
                            # append to buffer
                            if rname not in output[dir]:
                                output[dir][rname] = []
                            output[dir][rname].append(str(lines[dir][l][a]))
            return output

class Polyalign:
    """
    Object to handle separate exhaustive bwa alignments of left and right reads.
    Estimates insert size and returns sam of only accepted best pairings.
    """
    def __init__(self, subject_path, reads_1_path, reads_2_path, output_path=".", output_basename="polyalign", output_type="filtered", retain_unmapped=None, read_length=None, read_orientation=None, insert_size_range=None, bwa_options={}, do_bwa_indexing=True, bwa_workers=None, python_workers=None):
        self.subject_path = subject_path
        self.reads_1_path = reads_1_path
        self.reads_2_path = reads_2_path
        self.output_path = output_path
        self.output_basename = output_basename
        self.output_type = output_type
        # set workers
        self.bwa_workers = bwa_workers
        self.python_workers = python_workers
        # set necessary bwa options, -a and -Y flags
        bwa_options["a"] = None
        bwa_options["Y"] = None
        # initialise bwa mem
        self.bwa_mem = BwaMem(subject_path, bwa_options, do_indexing=do_bwa_indexing, cpus=bwa_workers)
        # set read orientation and insert size, fetching automatically if necessary
        self.read_length = read_length
        self.read_orientation = read_orientation
        self.insert_size_range = insert_size_range
        if self.read_length is None or self.read_orientation is None or self.insert_size_range is None:
            self.estimateInsertSizeRange()
        self.retain_unmapped = retain_unmapped
        if self.retain_unmapped is None and output_type == "paired":
            self.retain_unmapped = False
        else:
            self.retain_unmapped = True

    def estimateInsertSizeRange(self, read_sample = 100000, insert_size_percentile = 0.01):
        """
        Using a sample of read_sample reads, estimate insert size distribution.
        Return an acceptable size range of 1st to 99th percentile.
        """
        def sampleReads(input_path, output_path, read_sample):
            with open(input_path, "r") as input_file:
                with open(output_path, "w") as output_file:
                    count = 0
                    while count < read_sample * 4:
                        line = input_file.readline()
                        if not line:
                            break
                        output_file.write(line)
                        count += 1
        
        print("[PA::EISR]", "Getting mate pair properties from read subset alignment")
        # output stats
        insert_size = []
        orientation = {"fr": 0, "ff": 0, "rr": 0, "rf": 0}
        read_length = []
        # sample reads
        print("[PA::EISR]", "Sampling "+str(read_sample)+" reads")
        sampleReads(self.reads_1_path, "tmp1.fastq.tmp", read_sample)
        sampleReads(self.reads_2_path, "tmp2.fastq.tmp", read_sample)
        # do alignment
        sam_1 = self.bwa_mem.mem("tmp1.fastq.tmp")
        sam_2 = self.bwa_mem.mem("tmp2.fastq.tmp")
        # get stats
        # read and ignore header lines
        sam_1.readLine()
        sam_1.readHeader()
        sam_2.readLine()
        sam_2.readHeader()
        # iterate through sets of alignment for each read
        next_read_1 = sam_1.nextRead()
        next_read_2 = sam_2.nextRead()
        PairFilter = AlignmentPairFilter()
        while next_read_1 is not None and next_read_2 is not None:
            next_read_1 = [Alignment(line) for line in next_read_1]
            next_read_2 = [Alignment(line) for line in next_read_2]
            # check for read name matching errors
            if next_read_1[0].qname != next_read_2[0].qname:
                print("[PA::EISR]", "Error: Read names do not match")
                raise ValueError("Read names do not match, are the reads the same order in the two FASTQ files?")
            # if a single pair
            if PairFilter.isSinglePair(next_read_1, next_read_2):
                # get and record stats
                read_pair = AlignmentPair(next_read_1[0], next_read_2[0])
                if read_pair.sameReference():
                    insert_size.append(read_pair.insertSize())
                    orientation[read_pair.orientation()] += 1
                    read_length += read_pair.readLength()
            next_read_1 = sam_1.nextRead()
            next_read_2 = sam_2.nextRead()
        # analyse insert_size
        insert_size.sort()
        insert_lower = insert_size[int(len(insert_size) * insert_size_percentile)]
        insert_upper = insert_size[int(len(insert_size) * (1 - insert_size_percentile))]
        read_length = sum(read_length) // len(read_length)
        if len(insert_size) == 0:
            print("[PA::EISR]", "Error: No valid pairs found from " + str(read_sample) + " sampled reads")
            raise ValueError("No valid pairs found")
        if len(insert_size) < read_sample / 2:
            print("[PA::EISR]", "Warning: Fewer than 50% of reads have valid pairs")
        if insert_lower >= insert_upper or insert_lower < read_length:
            print("[PA::EISR]", "Error: Insert size range is invalid")
            print("[PA::EISR]", "Read length:", read_length)
            print("[PA::EISR]", "Insert size range:", insert_lower, insert_upper)
            raise ValueError("Insert size range is invalid")
        # print result
        print("[PA::EISR]", "Result:")
        print("[PA::EISR]", "Read length:", read_length)
        print("[PA::EISR]", "Insert size "+str(insert_size_percentile * 100)+"% and "+str(100 - insert_size_percentile * 100)+"% percentile:", insert_lower, insert_upper)
        print("[PA::EISR]", "Orientation:", ", ".join([x[0]+": "+str(x[1]) for x in orientation.items()]))
        # set properties
        if self.read_orientation is None:
            self.read_orientation = max(orientation, key=orientation.get)
        if self.insert_size_range is None:
            self.insert_size_range = (insert_lower, insert_upper)
        if self.read_length is None:
            self.read_length = read_length
        # remove temporary files
        os.remove("tmp1.fastq.tmp")
        os.remove("tmp2.fastq.tmp")

    def batchChunkWorker(self, chunk):
        """
        Run pairing analysis on list of reads
        """
        start_time = time.time()
        # parse lines as alignments
        chunk["list_1"] = [[Alignment(line) for line in read] for read in chunk["list_1"]]
        chunk["list_2"] = [[Alignment(line) for line in read] for read in chunk["list_2"]]
        # do filtering
        PairFilter = AlignmentPairFilter(chunk["read_length"], chunk["read_orientation"], chunk["insert_size_range"], chunk["retain_unmapped"], chunk["output_type"])
        current_stats = {}
        result_1, result_2 = [], []
        for index in range(len(chunk["list_1"])):
            curr_result_1, curr_result_2, status = PairFilter.analyseReadPairing(chunk["list_1"][index], chunk["list_2"][index])
            if curr_result_1 is not None and curr_result_2 is not None and status is not None:
                result_1.append(curr_result_1)
                result_2.append(curr_result_2)
                if status not in current_stats:
                    current_stats[status] = 0
                current_stats[status] += 1
        output = PairFilter.alignmentOutputPreprocess(result_1, result_2)
        # return
        chunk["time"]["process time"] = time.time() - start_time
        return {
            "output": output,
            "index": chunk["index"],
            "time": chunk["time"],
            "stats": current_stats,
        }

    def polyalign(self, mode="parallel", chunk_length=1000):
        """
        Do the polyalignment.
        Aligns left and right reads separately, parsing output for all possible correct pairings.
        """
        def returnHeader():
            """
            Read input sams and return header lines
            """
            command = sys.argv.copy()
            command[0] = os.path.basename(command[0])
            program_line = Header(["@PG", "ID:polyalign", "PN:polyalign", "VN:0.0.0", "CL:"+" ".join(command)])
            return [Header(line) for line in sam_1.readHeader()] + [program_line], [Header(line) for line in sam_2.readHeader()] + [program_line]

        print("[PA::PA]", "Polyaligning in", self.output_type, "mode")
        # set mode
        if mode is None:
            mode = "parallel"
        if mode not in ["serial", "parallel"]:
            raise ValueError(f"Unknown mode '{mode}'")
        print("[PA::PA]", "Parallelisation mode:", mode)
        # initialise stats
        self.stats = {}
        read_index = 0
        # open output
        output = Output(self.output_basename, self.output_type, self.output_path)
        # start alignment
        sam_1 = self.bwa_mem.mem(self.reads_1_path)
        sam_2 = self.bwa_mem.mem(self.reads_2_path)
        # initialise input
        sam_1.readLine()
        sam_2.readLine()
        # header
        header_1, header_2 = returnHeader()
        output.writeHeader(header_1, header_1)
        # alignment and filtering
        # set multiprocess parameters
        chunk_reporting_interval = max(1, 100000//chunk_length)
        multiprocess_mode = "process"
        workers = self.python_workers
        if workers is None:
            workers = multiprocessing.cpu_count()
        else:
            workers = min(max(1, workers), multiprocessing.cpu_count())
        if mode == "serial" or workers == 1:
            multiprocess_mode = None
            workers = 1
        # setup executor for multiprocessing mode
        if multiprocess_mode is None:
            # run in a single thread, single ThreadPoolExecutor is equivalent to simple ongoing script
            Executor = concurrent.futures.ThreadPoolExecutor
            workers = 1
        elif multiprocess_mode == "process":
            Executor = concurrent.futures.ProcessPoolExecutor
        elif multiprocess_mode == "thread":
            Executor = concurrent.futures.ThreadPoolExecutor
        else:
            raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")
        print("[PA::RBA]", "Workers:", workers, ", Chunk length:", chunk_length, ", Multiprocessing mode:", multiprocess_mode)
        # start execution
        # helper function for prepping a chunk
        def prep_chunk(index, sam_1, sam_2, chunk_length):
            start_time = time.time()
            chunk = {
                "index": index,
                "time": {},
                "list_1": [],
                "list_2": [],
                "read_length": self.read_length,
                "read_orientation": self.read_orientation,
                "insert_size_range": self.insert_size_range,
                "retain_unmapped": self.retain_unmapped,
                "output_type": self.output_type,
            }
            for _ in range(chunk_length):
                next_read_1 = sam_1.nextRead()
                next_read_2 = sam_2.nextRead()
                if next_read_1 is not None and next_read_2 is not None:
                    chunk["list_1"].append(next_read_1)
                    chunk["list_2"].append(next_read_2)
            chunk["time"]["read time"] = time.time() - start_time
            if len(chunk["list_1"]) > 0 and len(chunk["list_2"]) > 0:
                return chunk
            else:
                print("[PA::PA]", "End of bwa mem output reached")
                return None
        # helper function for handling output
        def handle_output(result):
            start_time = time.time()
            for status in result["stats"]:
                if status not in self.stats:
                    self.stats[status] = 0
                self.stats[status] += result["stats"][status]
            output.writePreprocessed(result["output"])
            result["time"]["write time"] = time.time() - start_time
            if result["index"]%chunk_reporting_interval == 0:
                print("[PA::PA]", "Cumulative stats at chunk", str(result["index"])+":", ", ".join([x[0]+": "+str(x[1]) for x in self.stats.items()]))
                print("[PA::PA]", "Performance on chunk", str(result["index"])+":", ", ".join([x[0]+": "+str(round(x[1], 5))+"s" for x in result["time"].items()]))
        # do analysis
        chunk_index = 0
        if workers > 1:
            # fire up executor
            with Executor(workers) as executor:
                concurrent_chunks = workers * 2
                futures = []
                # start initial chunks, up to concurrent_chunks
                chunk = prep_chunk(chunk_index, sam_1, sam_2, chunk_length)
                while chunk is not None and chunk_index < concurrent_chunks:
                    futures.append(executor.submit(self.batchChunkWorker, chunk=chunk))
                    chunk_index += 1
                    chunk = prep_chunk(chunk_index, sam_1, sam_2, chunk_length)
                # as long as chunks can be assembled, start new chunks when tasks complete
                while chunk is not None:
                    # check which chunks are complete
                    complete = []
                    for future_index in range(len(futures)-1, -1, -1):
                        if futures[future_index].done():
                            complete.append(future_index)
                    # submit that many new chunks
                    i = 0
                    while chunk is not None and i < len(complete):
                        futures.append(executor.submit(self.batchChunkWorker, chunk=chunk))
                        chunk_index += 1
                        chunk = prep_chunk(chunk_index, sam_1, sam_2, chunk_length)
                        i += 1
                    # handle output from completed chunks
                    for future_index in complete:
                        future = futures.pop(future_index)
                        handle_output(future.result())
                    time.sleep(0.001)
                # wait for all remaining futures to complete and handle output
                for future in concurrent.futures.as_completed(futures):
                    handle_output(future.result())
        else:
            # single thread
            chunk = prep_chunk(0, sam_1, sam_2, chunk_length)
            while chunk is not None:
                handle_output(self.batchChunkWorker(chunk=chunk))
                chunk_index += 1
                chunk = prep_chunk(chunk_index, sam_1, sam_2, chunk_length)
        # end output
        output.closeOutput()
        print("[PA::PA]", "Stats:", ", ".join([x[0]+": "+str(x[1]) for x in self.stats.items()]))
        print("[PA::PA]", "Polyaligning complete")

class Output:
    def __init__(self, base_name, output_mode, base_path, append_output=False):
        self.base_path = base_path
        self.base_name = base_name
        self.output_mode = output_mode
        self.append_output = append_output # whether to remove existing output (False) files or happily append (True)
        self.soft_file_limit = None
        output_modes = ["filtered", "paired", "splitfiltered"]
        if self.output_mode not in output_modes:
            raise ValueError("Output mode must be one of:", ", ".join(output_modes))
        if self.base_name == "-" and output_mode != "filtered":
            raise ValueError("Output mode except for 'paired' requires a base path")
        self.output_paths = self.setupOutput()
        print("[PA::OP]", "Output mode:", self.output_mode, "Base name:", self.base_name, "Base path:", self.base_path)

    def setupOutput(self):
        # setup output files
        mode = "w"
        if self.append_output:
            mode = "a"
        if self.base_name == "-":
            print("[PA::OP]", "Outputting to stdout")
            return [sys.stdout]
        if self.output_mode == "paired":
            print("[PA::OP]", "Output file:", self.base_name+".sam")
            if not self.append_output:
                if os.path.exists(os.path.join(self.base_path, self.base_name+".sam")):
                    os.remove(os.path.join(self.base_path, self.base_name+".sam"))
            return [open(os.path.join(self.base_path, self.base_name + ".sam"), mode)]
        if self.output_mode == "filtered":
            print("[PA::OP]", "Output files:", self.base_name+"_1.sam", self.base_name+"_2.sam")
            if not self.append_output:
                if os.path.exists(os.path.join(self.base_path, self.base_name+"_1.sam")):
                    os.remove(os.path.join(self.base_path, self.base_name+"_1.sam"))
                if os.path.exists(os.path.join(self.base_path, self.base_name+"_2.sam")):
                    os.remove(os.path.join(self.base_path, self.base_name+"_2.sam"))
            return [open(os.path.join(self.base_path, self.base_name + "_1.sam"), mode), open(os.path.join(self.base_path, self.base_name + "_2.sam"), mode)]
        if self.output_mode == "splitfiltered":
            print("[PA::OP]", "Output directories:", self.base_name+"_1", self.base_name+"_2")
            if not self.append_output:
                self.checkOutputDir(os.path.join(self.base_path, self.base_name+"_1"))
                self.checkOutputDir(os.path.join(self.base_path, self.base_name+"_2"))
            self.soft_file_limit = self.getOpenFileLimit()
            self.open_files = [{}, {}] # open file dicts, index 0 for read 1, 1 for read 2
            self.line_buffers = [{}, {}] # line buffers for each reference sequence, index 0 for read 1, 1 for read 2
            self.buffer_lines = 1000
            return [os.path.join(self.base_path, self.base_name+"_1"), os.path.join(self.base_path, self.base_name+"_2")]

    def getOpenFileLimit(self):
        """
        Get the maximum number of open files
        """
        softfilelimit = 1024
        import importlib.util
        if importlib.util.find_spec("resource") is not None:
            # if resource, a linux-only module, is available
            import resource
            softfilelimit, hardfilelimit = resource.getrlimit(resource.RLIMIT_NOFILE)
        print("[PA::OP]", "Soft file limit:", softfilelimit)
        return softfilelimit - 10

    def checkOutputDir(self, path):
        """
        Make the output path or empty it if it is not empty
        """
        # check the output directory and make if needed
        if os.path.exists(path):
            if len(os.listdir(path)) > 0:
                # remove all existing files
                for file in os.listdir(path):
                    os.remove(os.path.join(path, file))
        else:
            os.mkdir(path)
        
    def writeLine(self, file, line):
        if line is not None:
            if self.base_path == "-":
                print(str(line))
            else:
                file.write(str(line)+"\n")
    
    def writeLines(self, file, lines):
        if lines is not None:
            for line in lines:
                self.writeLine(file, line)
    
    def writeHeader(self, header_1, header_2):
        print("[PA::OP]", "Writing headers")
        if self.output_mode == "filtered":
            # two output files, so write header to both
            self.writeLines(self.output_paths[0], header_1)
            self.writeLines(self.output_paths[1], header_2)
        if self.output_mode == "paired":
            # single output file, so single header
            self.writeLines(self.output_paths[0], header_1)
        if self.output_mode == "splitfiltered":
            # add headers to line buffers for each reference sequence
            headers = [header_1, header_2]
            for h in range(len(headers)):
                for line in headers[h]:
                    if line.line[0] == "@SQ":
                        # for sequence lines, assumes they come first
                        for entry in line.line[1:]:
                            if entry[0:2] == "SN":
                                referencename = entry[3:]
                        # setup output buffer for each reference sequence, with @SQ specific to output file
                        self.line_buffers[h][referencename] = []
                        self.line_buffers[h][referencename].append(line)
                    else:
                        # all other header lines to all files
                        for reference in self.line_buffers[h]:
                            self.line_buffers[h][reference].append(line)

    def writePreprocessed(self, output):
        """
        Writes lines preprocessed by the alignmentOutputPreprocess function
        """
        def bufferedLineWrite(dir, rname, lines):
            # setup line buffer, if not already done for this rname
            if rname not in self.line_buffers[dir]:
                self.line_buffers[dir][rname] = []
            # add to line buffer
            self.line_buffers[dir][rname] += lines
            # write buffer to file if too many lines buffered
            if len(self.line_buffers[dir][rname]) > self.buffer_lines:
                # check if too many files are open, and close first opened one if necessary
                if len(self.open_files[dir]) >= self.soft_file_limit / 2:
                    # get the first open file and write any buffered lines for that rname to it
                    first_file = next(iter(self.open_files[dir]))
                    self.writeLines(self.open_files[dir][first_file], self.line_buffers[dir][first_file])
                    # close the first open file and remove from open files list
                    self.open_files[dir][first_file].close()
                    self.open_files[dir].pop(first_file)
                # open the output file if necessary
                if rname not in self.open_files[dir]:
                    self.open_files[dir][rname] = open(os.path.join(self.output_paths[dir], rname+"_"+str(dir+1)+".sam"), "a")
                # write lines and clear buffer
                self.writeLines(self.open_files[dir][rname], self.line_buffers[dir][rname])
                self.line_buffers[dir][rname] = []

        if self.output_mode == "filtered":
            self.writeLines(self.output_paths[0], output[0])
            self.writeLines(self.output_paths[1], output[1])
        if self.output_mode == "paired":
            self.writeLines(self.output_paths[0], output)
        if self.output_mode == "splitfiltered" or self.output_mode == "split":
            for rname in output[0]:
                bufferedLineWrite(0, rname, output[0][rname])
            for rname in output[1]:
                bufferedLineWrite(1, rname, output[1][rname])

    def closeOutput(self):
        if self.base_path == "-":
            return
        if self.output_mode == "filtered":
            self.output_paths[0].close()
            self.output_paths[1].close()
        if self.output_mode == "paired":
            self.output_paths[0].close()
        if self.output_mode == "splitfiltered":
            # write buffers for open files
            for l in range(len(self.line_buffers)):
                for file in self.line_buffers[l]:
                    if self.open_files[l][file] is not None:
                        self.writeLines(self.open_files[l][file], self.line_buffers[l][file])
                        self.line_buffers[l][file] = []
                        self.open_files[l][file].close()
            # write remaining buffers
            for l in range(len(self.line_buffers)):
                for reference in self.line_buffers[l]:
                    if len(self.line_buffers[l][reference]) > 0:
                        if self.open_files[l][reference] is None:
                            self.open_files[l][reference] = open(os.path.join(self.base_path, reference+"_"+str(l+1)+".sam"), "w")
                        self.writeLines(self.open_files[l][reference], self.line_buffers[l][reference])
                        self.line_buffers[l][reference] = []
                        self.open_files[l][reference].close()

class Splitfasta:
    """
    Object to split an input file into one fasta file per contig.
    Helper tool to prepare split fasta files, for use with splitfiltered alignments.
    """
    def __init__(self, fasta_path, output_path=".", output_basename="polyalign"):
        self.fasta_path = fasta_path
        self.output_path = output_path
        self.output_basename = output_basename
    
    def splitfasta(self):
        """
        Split fasta file into one file per contig
        """
        print("[PA::SF]", "Splitting fasta file into one file per contig")
        # make and check output directory
        if not os.path.exists(os.path.join(self.output_path, self.output_basename)):
            os.mkdir(os.path.join(self.output_path, self.output_basename))
        if len(os.listdir(os.path.join(self.output_path, self.output_basename))) > 0:
            raise ValueError("Output directory must be empty")
        # open input file
        with open(self.fasta_path, "r") as input_file:
            contig_name = None
            output_file = None
            # read first line, then iterate
            line = input_file.readline()
            while line:
                # if line is a header, setup new output and write header
                if line[0] == ">":
                    if output_file is not None:
                        output_file.close()
                    contig_name = line[1:].strip().split()[0]
                    output_file = open(os.path.join(self.output_path, self.output_basename, contig_name+".fasta"), "w")
                    output_file.write(line)
                else:
                    output_file.write(line)
                line = input_file.readline()
            # close final output file
            if output_file is not None:
                output_file.close()
        print("[PA::SF]", "Splitting complete")
