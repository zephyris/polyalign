import os, sys, subprocess, random, multiprocessing, concurrent.futures, copy
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
            self.cpus = multiprocessing.cpu_count() / 2
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
            print("-PA-", "Error: " + self.program + " not found")

    def makeIndex(self):
        """
        Index subject_path fasta file using bwa index.
        """
        index_exists = True
        for suffix in self.index_suffices:
            if not os.path.exists(self.subject_path + suffix):
                index_exists = False
                print("-PA-", "Index files not found,", self.subject_path + suffix, "missing")
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
            return [Header(line) for line in lines]
        
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
            return [Alignment(line) for line in lines]
    
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
        print("-PA-", "bwa mem command:", " ".join(command))
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
        self.qname = sam_line[0]
        self.flag = int(sam_line[1])
        self.rname = sam_line[2]
        self.pos = int(sam_line[3])
        self.mapq = int(sam_line[4])
        self.cigar = sam_line[5]
        self.rnext = sam_line[6]
        self.pnext = int(sam_line[7])
        self.tlen = int(sam_line[8])
        self.seq = sam_line[9]
        self.qual = sam_line[10]
        self.optional = sam_line[11:]
        self.flag_indices = ["primary", "properly_aligned", "unmapped", "mate_unmapped", "reverse", "mate_reverse", "first", "second", "secondary", "qc_fail", "duplicate", "supplementary"]
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """
        Return sam line as a string
        """
        return "\t".join(self.columns())
    
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
        start, end = self.pos - 1, self.pos - 1 + self.tlen
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

class Polyalign:
    """
    Object to handle separate exhaustive bwa alignments of left and right reads.
    Estimates insert size and returns sam of only accepted best pairings.
    """
    def __init__(self, subject_path, reads_1_path, reads_2_path, output_path=".", output_basename="polyalign", output_type="filtered", read_orientation=None, insert_size_range=None, bwa_options={}, do_bwa_indexing=True):
        self.subject_path = subject_path
        self.reads_1_path = reads_1_path
        self.reads_2_path = reads_2_path
        self.output_path = output_path
        self.output_basename = output_basename
        output_types = ["filtered", "paired"]
        self.output_type = output_type
        if self.output_type not in output_types:
            raise ValueError("Output type must be one of: "+", ".join(output_types))
        # set necessary bwa options, -a and -Y flags
        bwa_options["a"] = None
        bwa_options["Y"] = None
        # initialise bwa mem
        self.bwa_mem = BwaMem(subject_path, bwa_options, do_indexing=do_bwa_indexing)
        # set read orientation and insert size, fetching automatically if necessary
        self.read_orientation = read_orientation
        self.insert_size_range = insert_size_range
        if self.read_orientation is None or self.insert_size_range is None:
            self.estimateInsertSizeRange()

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
    
    def findGoodPairs(self, read_1, read_2):
        """
        Find all good pairs of read_1 and read_2 alignments, pairs which align to the same reference, with the correct insert size, in the correct orientation
        """
        good_1, good_2 = [False] * len(read_1), [False] * len(read_2)
        paired_1, paired_2 = [], []
        for index_1, alignment_1 in enumerate(read_1):
            for index_2, alignment_2 in enumerate(read_2):
                if self.isGoodSinglePair([alignment_1], [alignment_2]):
                    good_1[index_1] = True
                    good_2[index_2] = True
                    paired_1.append(alignment_1), paired_2.append(alignment_2)
        if self.output_type == "filtered":
            return [read_1[index] for index, value in enumerate(good_1) if value], [read_2[index] for index, value in enumerate(good_2) if value]
        if self.output_type == "paired":
            return paired_1, paired_2
    
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
        
        print("-PA-", "Getting mate pair properties from read subset alignment")
        # output stats
        insert_size = []
        orientation = {"fr": 0, "ff": 0, "rr": 0, "rf": 0}
        # sample reads
        print("-PA-", "Sampling "+str(read_sample)+" reads")
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
        while next_read_1 is not None:
            # if a single pair
            if self.isSinglePair(next_read_1, next_read_2):
                # get and record statsenumerate
                read_pair = AlignmentPair(next_read_1[0], next_read_2[0])
                if read_pair.sameReference():
                    insert_size.append(read_pair.insertSize())
                    orientation[read_pair.orientation()] += 1
            next_read_1 = sam_1.nextRead()
            next_read_2 = sam_2.nextRead()
        # analyse insert_size
        insert_size.sort()
        insert_lower = insert_size[int(len(insert_size) * insert_size_percentile)]
        insert_upper = insert_size[int(len(insert_size) * (1 - insert_size_percentile))]
        # print result
        print("-PA-", "Result:")
        print("-PA-", "Insert size "+str(insert_size_percentile * 100)+"% and "+str(100 - insert_size_percentile * 100)+"% percentile:", insert_lower, insert_upper)
        print("-PA-", "Orientation:", ", ".join([x[0]+": "+str(x[1]) for x in orientation.items()]))
        # set properties
        if self.read_orientation is None:
            self.read_orientation = max(orientation, key=orientation.get)
        if self.insert_size_range is None:
            self.insert_size_range = (insert_lower, insert_upper)
        # remove temporary files
        os.remove("tmp1.fastq.tmp")
        os.remove("tmp2.fastq.tmp")

    def polyalign(self, retain_unmapped=False, mode="parallel"):
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
            return sam_1.readHeader() + [program_line], sam_2.readHeader() + [program_line]
        
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
                print("-PA-", "Error: Read names do not match")
            # check for no alignments
            if self.isUnalignedPair(next_read_1, next_read_2):
                # For pairs where both have zero alignments, retain sam lines, but no way to contribute to polishing.
                status = "both unaligned"
                if retain_unmapped:
                    return next_read_1, next_read_2, status
                return [], [], status
            if self.isOneUnalignedPair(next_read_1, next_read_2):
                # For pairs with exactly one alignment of each read, keep sam lines, ie. treat unique paired alignments as correct.
                status = "single unaligned"
                return next_read_1, next_read_2, status
            if self.isSinglePair(next_read_1, next_read_2):
                # For pairs where one has zero alignments, keep sam lines, ie. treat single alignments as correct.
                status = "unique both aligned"
                return next_read_1, next_read_2, status
            # For pairs with multiple alignments, keep sam lines if it pairs (same ref seq, correct ori, good insert) with any of its pair's alignments.
            good_1, good_2 = self.findGoodPairs(next_read_1, next_read_2)
            if len(good_1) == 0 or len(good_2) == 0:
                status = "no good pairs from multiple both aligned"
                return [], [], status
            # add sequence and quality to first good reads, replace with "*" for the rest
            for index, good in enumerate(good_1):
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
            status = "good pairs from multiple both aligned"
            return good_1, good_2, status

        def batchChunkWorker(chunk=None):
            """
            Run analysis_function on list of reads, using multiprocessing or multithreading
            """
            current_stats = {}
            result_1, result_2 = [], []
            for index in range(len(chunk["list_1"])):
                curr_result_1, curr_result_2, status = analyseReadPairing(self, chunk["list_1"][index], chunk["list_2"][index])
                if curr_result_1 is not None and curr_result_2 is not None and status is not None:
                    result_1 += curr_result_1
                    result_2 += curr_result_2
                    if status not in current_stats:
                        current_stats[status] = 0
                    current_stats[status] += 1
            result = {
                "index": chunk["index"],
                "chunk": (result_1, result_2),
                "stats": current_stats
            }
            return result

        def readBatchAnalysis(chunk_length=100000, multiprocess_mode="thread", workers=None):
            """
            Run analysis_function on list of reads, return processed list
            """
            # set number of workers
            if workers is None:
                workers = multiprocessing.cpu_count()
            # setup multiprocessing mode
            if multiprocess_mode is None:
                # run in a single thread, still use the ThreadPoolExecutor since that's equivalent
                print("-PA-", "Single thread")
                Executor = concurrent.futures.ThreadPoolExecutor
                workers = 1
            elif multiprocess_mode == "process":
                print("-PA-", "Parallel processes with", workers, "workers")
                # setup executor as a process pool
                Executor = concurrent.futures.ProcessPoolExecutor
            elif multiprocess_mode == "thread":
                # setup executor as a thread pool
                print("-PA-", "Parallel threads with", workers, "workers")
                Executor = concurrent.futures.ThreadPoolExecutor
            else:
                raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")
            # fetch data for worklist
            list_chunks = []
            for i in range(workers):
                list_chunks.append({"index": i, "list_1": [], "list_2": []})
                for j in range(chunk_length):
                    next_read_1 = sam_1.nextRead()
                    next_read_2 = sam_2.nextRead()
                    if next_read_1 is not None and next_read_2 is not None:
                        list_chunks[i]["list_1"].append(next_read_1)
                        list_chunks[i]["list_2"].append(next_read_2)
            print("-PA-", "Analysing read pairing,", len(list_chunks), "chunks of", len(list_chunks[0]["list_1"]), "read pairs")
            # fire up executor
            with Executor(workers) as executor:
                futures = [executor.submit(batchChunkWorker, chunk=chunk) for chunk in list_chunks]
                results = [future.result() for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), smoothing=0, disable=True)]
            # reconstruct stats
            for result in results:
                for key in result["stats"]:
                    if key not in self.stats:
                        self.stats[key] = 0
                    self.stats[key] += result["stats"][key]
            # sort by chunk index and reconstruct results
            results.sort(key=lambda x: x["index"])
            result_1, result_2 = [], []
            for result in results:
                result_1 += result["chunk"][0]
                result_2 += result["chunk"][1]
            return (result_1, result_2)
        
        def printLine(line):
            if line is not None:
                print(str(line))
        
        def printLines(lines):
            if lines is not None:
                for line in lines:
                    printLine(line)
        
        def writeLine(file, line):
            if line is not None:
                file.write(str(line)+"\n")
        
        def writeLines(file, lines):
            if lines is not None:
                for line in lines:
                    writeLine(file, line)

        print("-PA-", "Polyaligning in", self.output_type, "mode")
        # set mode
        mode = "parallel"
        if mode is None:
            mode = "parallel"
        if mode not in ["serial", "parallel"]:
            raise ValueError(f"Unknown mode '{mode}'")
        # initialise stats
        self.stats = {}
        read_index = 0
        # open output
        if self.output_type == "filtered":
            file_1 = open(os.path.join(self.output_path, self.output_basename+"_1.sam"), "w")
            file_2 = open(os.path.join(self.output_path, self.output_basename+"_2.sam"), "w")
            print("-PA-", "Output files:", self.output_basename+"_1.sam", self.output_basename+"_2.sam")
        elif self.output_type == "paired":
            if self.output_basename == "-":
                file = sys.stdout
                print("-PA-", "Outputting to stdout")
            else:
                file = open(os.path.join(self.output_path, self.output_basename+".sam"), "w")
                print("-PA-", "Output file:", self.output_basename+".sam")
        # start alignment
        sam_1 = self.bwa_mem.mem(self.reads_1_path)
        sam_2 = self.bwa_mem.mem(self.reads_2_path)
        # initialise input
        sam_1.readLine()
        sam_2.readLine()
        # header
        header_1, header_2 = returnHeader()
        if self.output_type == "filtered":
            writeLines(file_1, header_1)
            writeLines(file_2, header_2)
        elif self.output_type == "paired":
            writeLines(file, header_1)
        # alignment and filtering
        if mode == "parallel":
            results = readBatchAnalysis()
            while len(results[0]) > 0 and len(results[1]) > 0:
                print("-PA-", "Stats:", ", ".join([x[0]+": "+str(x[1]) for x in self.stats.items()]))
                if self.output_type == "filtered":
                    writeLines(file_1, results[0])
                    writeLines(file_2, results[1])
                elif self.output_type == "paired":
                    # TODO: Output is currently simple interleaving, but should set proper values...
                    # set [6] (rnext), [7] (pnext), sam flags [1]
                    writeLines(file, [x for a in zip(results[0], results[1]) for x in a])
                results = readBatchAnalysis()
        elif mode == "serial":
            next_read_1 = sam_1.nextRead()
            next_read_2 = sam_2.nextRead()
            alignment_1, alignment_2, status = analyseReadPairing(self, next_read_1, next_read_2)
            while alignment_1 is not None and alignment_2 is not None:
                read_index += 1
                if status is not None:
                    if status not in self.stats:
                        self.stats[status] = 0
                    self.stats[status] += 1
                if read_index % 100000 == 0:
                    print("-PA-", "Stats:", ", ".join([x[0]+": "+str(x[1]) for x in self.stats.items()]))
                if self.output_type == "filtered":
                    writeLines(file_1, alignment_1)
                    writeLines(file_2, alignment_2)
                elif self.output_type == "paired":
                    # TODO: Output is currently simple interleaving, but should set proper values...
                    # set [6] (rnext), [7] (pnext), sam flags [1]
                    writeLines(file, [x for a in zip(alignment_1, alignment_2) for x in a])
                next_read_1 = sam_1.nextRead()
                next_read_2 = sam_2.nextRead()
                alignment_1, alignment_2, status = analyseReadPairing(self, next_read_1, next_read_2)
        # end output
        print("-PA-", "End of bwa mem output reached")
        print("-PA-", "Final stats:", ", ".join([x[0]+": "+str(x[1]) for x in self.stats.items()]))
        print("-PA-", "Polyaligning complete")
        if self.output_type == "filtered":
            file_1.close()
            file_2.close()
        elif self.output_type == "paired":
            if self.output_basename != "-":
                file.close()
