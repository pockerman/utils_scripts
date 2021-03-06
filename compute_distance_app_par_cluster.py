import itertools
import csv
import multiprocessing as mp
import json
import sys
import textdistance

INFO = "INFO: "


class NuclOutFileReader(object):

    def __init__(self, exclude_seqs=["NO_REPEATS"], delimiter: str='\t') -> None:
        self._delimiter = delimiter
        self._exclude_seqs = exclude_seqs

    def read_line_nucl_out_file(self, line: str, delimiter='\t') -> tuple:
        line_data = line.split(delimiter)
        chromosome = line_data[0].strip()
        start = int(line_data[1])
        end = int(line_data[2])
        seq = line_data[3].strip()
        state = line_data[4].strip()
        return chromosome, start, end, seq, state

    def __call__(self, filename) -> list:
        with open(filename, 'r', newline="\n") as fh:
            seqs = []

            for line in fh:
                chromosome, start, end, seq, state = self.read_line_nucl_out_file(line=line, delimiter=self._delimiter)

                if seq not in self._exclude_seqs:
                    seqs.append([chromosome, start, end, seq, state])

            return seqs


class TextDistanceCalculator(object):
    """
    Wrapper class for text distance calculation
    """

    NAMES = ['ham', 'mlipns', 'lev', 'damlev', 'jwink', 'str',
             'nw', 'sw', 'got', 'jac', 'sor', 'tve', 'ov', 'tan',
             'cos', 'mon', 'bag', 'lcsseq', 'lcsstr', 'rat',
             'ari', 'rle', 'bwt', 'sqr', 'ent', 'bz2', 'lzm',
             'zli', 'mra', 'edi', 'pre', 'pos',
             'len', 'id', 'mat',  'CPF', 'all']

    @staticmethod
    def build_calculator(name):

        if name not in TextDistanceCalculator.NAMES:
            raise Error("Distance type '{0}' is invalid".format(name))

        if name   == 'ham'   : return textdistance.Hamming()
        elif name == 'mlipns': return textdistance.MLIPNS()
        elif name == 'lev'   : return textdistance.Levenshtein()
        elif name == 'damlev': return textdistance.DamerauLevenshtein()
        elif name == 'jwink' : return textdistance.JaroWinkler()
        elif name == 'str'   : return textdistance.StrCmp95()
        elif name == 'nw'    : return textdistance.NeedlemanWunsch()
        elif name == 'sw'    : return textdistance.SmithWaterman()
        elif name == 'got'   : return textdistance.Gotoh()
        elif name == 'jac'   : return textdistance.Jaccard()
        elif name == 'sor'   : return textdistance.Sorensen()
        elif name == 'tve'   : return textdistance.Tversky()
        elif name == 'ov'    : return textdistance.Overlap()
        elif name == 'tan'   : return textdistance.Tanimoto()
        elif name == 'cos'   : return textdistance.Cosine()
        elif name == 'mon'   : return textdistance.MongeElkan()
        elif name == 'bag'   : return textdistance.Bag()
        elif name == 'lcsseq': return textdistance.LCSSeq()
        elif name == 'lcsstr': return textdistance.LCSStr()
        elif name == 'rat'   : return textdistance.RatcliffObershelp()
        elif name == 'ari'   : return textdistance.ArithNCD()
        elif name == 'rle'   : return textdistance.RLENCD()
        elif name == 'bwt'   : return textdistance.BWTRLENCD()
        elif name == 'sqr'   : return textdistance.SqrtNCD()
        elif name == 'ent'   : return textdistance.EntropyNCD()
        elif name == 'bz2'   : return textdistance.BZ2NCD()
        elif name == 'lzm'   : return textdistance.LZMANCD()
        elif name == 'zli'   : return textdistance.ZLIBNCD()
        elif name == 'mra'   : return textdistance.MRA()
        elif name == 'edi'   : return textdistance.Editex()
        elif name == 'pre'   : return textdistance.Prefix()
        elif name == 'pos'   : return textdistance.Postfix()
        elif name == 'len'   : return textdistance.Length()
        elif name == 'id'    : return textdistance.Identity()
        elif name == 'mat'   : return textdistance.Matrix()
        elif name == "L2Norm": return L2Norm()
        elif name == 'CPF'   : return CPF()
        elif name == 'all'   : return AllDistancesCalculator()

    @staticmethod
    def read_sequence_comparison_file(filename, strip_path, delim=',', commment_delim='#'):

        with open(filename, 'r') as f:
            similarity_map = {}

            for line in f:

                if line.startswith(commment_delim):
                    continue

                line = line.split(delim)

                if strip_path:
                    line[0] = line[0].split('/')[-1]
                    line[1] = line[1].split('/')[-1]
                similarity_map[(line[0], line[1])] = float(line[2])

            return similarity_map

    def __init__(self, dist_type):

        if dist_type not in TextDistanceCalculator.NAMES:
            raise ValueError("Distance type '{0}' is invalid".format(dist_type))

        self._dist_type = dist_type

    def calculate(self, txt1, txt2, **options):

        # build a calculator
        calculator = TextDistanceCalculator.build_calculator(name=self._dist_type)

        set_options = getattr(calculator, "set_options", None)

        if set_options is not None:
            calculator.set_options(**options)

        return calculator.similarity(txt1, txt2)


def read_json(filename):

    """
    Read the json configuration file and
    return a map with the config entries
    """
    with open(filename) as json_file:
        json_input = json.load(json_file)
        return json_input


def reverse_complement_table(seq, tab):
    return seq.translate(tab)[::-1]


def write_pair_segments_distances(proc_id, master_proc, start, end, input_file,
                                  line_counter, outdir, distance_type, error_map):

    try:

        import sys
        from pathlib import Path

        # read original seqs
        filereader = NuclOutFileReader(exclude_seqs=["NO_REPEATS"])  # NuclOutFileReader()
        sequences = filereader(filename=input_file)

        if end == -1:
            end = len(sequences)

        if start >= end:
            raise ValueError("Error: start index={0} greater than end index={1}".format(start, end))

        if proc_id == master_proc:
            print("{0} Master process is working on [{1}, {2})".format(INFO, start, end))

        tab = str.maketrans("ACTGRMW", "TGACYKS")
        tabPP = str.maketrans("AGCT", "RRYY")
        tabAK = str.maketrans("ACGT", "MMKK")
        tabWS = str.maketrans("ATCG", "WWSS")

        # build the wrapper
        calculator = TextDistanceCalculator.build_calculator(name=distance_type)

        part_counter = 0
        touched = []
        lines = []

        # only work on the part assigned [start, end)
        for i in range(start, end, 1):
            for j in range(len(sequences)):
                if j >= start and (i, j) not in touched and (j, i) not in touched:

                    chr_seq_1 = sequences[i][0].strip()
                    start1 = sequences[i][1]
                    end1 = sequences[i][2]
                    seq1 = sequences[i][3].strip()
                    state1 = sequences[i][4].strip()

                    # print(sequences[j])
                    chr_seq_2 = sequences[j][0].strip()
                    start2 = sequences[j][1]
                    end2 = sequences[j][2]
                    seq2 = sequences[j][3].strip()
                    state2 = sequences[j][4].strip()

                    if len(seq1) + len(seq2) < 200:

                        # calculate the distance
                        distance1 = calculator.normalized_distance(seq1, seq2)
                        distance2 = calculator.normalized_distance(seq1,
                                                                   reverse_complement_table(seq=seq2, tab=tab))

                        distance = min(distance1, distance2)

                        PPseq1 = seq1.translate(tabPP)
                        AKseq1 = seq1.translate(tabAK)
                        WSseq1 = seq1.translate(tabWS)

                        PPseq2 = seq2.translate(tabPP)
                        AKseq2 = seq2.translate(tabAK)
                        WSseq2 = seq2.translate(tabWS)

                        PPdist = calculator.normalized_distance(PPseq1, PPseq2)
                        PPdist_rev = calculator.normalized_distance(PPseq1,
                                                                    reverse_complement_table(seq=PPseq2, tab=tab))

                        pp_distance = min(PPdist, PPdist_rev)

                        AKdist = calculator.normalized_distance(AKseq1, AKseq2)
                        AKdist_rev = calculator.normalized_distance(AKseq1,
                                                                    reverse_complement_table(seq=AKseq2, tab=tab))

                        ak_distance = min(AKdist, AKdist_rev)

                        WSdist = calculator.normalized_distance(WSseq1, WSseq2)
                        WSdist_rev = calculator.normalized_distance(WSseq1,
                                                                    reverse_complement_table(seq=WSseq2, tab=tab))

                        ws_distance = min(WSdist, WSdist_rev)

                        lines.append([chr_seq_1, start1, end1, seq1, state1,
                                      chr_seq_2, start2, end2, seq2, state2,
                                      distance, pp_distance, ak_distance, ws_distance])

                        # if we reached the batch then flush
                        # to the output file
                        if len(lines) == line_counter:

                            filename = Path(outdir + "part_seqs_pairs_proc_" + str(proc_id) + "_" + str(part_counter) + ".csv")

                            if filename.is_file():
                                raise ValueError("ERROR: filename={0} Exists".format(filename))

                            if proc_id == master_proc:
                                print("{0} Writing to {1}".format(INFO, filename))

                            with open(filename, 'w', newline="\n") as fh:
                                writer = csv.writer(fh, delimiter=",")
                                writer.writerow(["#", "ChrSeq-1", "StartSeq-1", "EndSeq-1", "Seq-1", "HMM-State-1",
                                                 "ChrSeq-2", "StartSeq-2", "EndSeq-2", "Seq-2", "HMM-State-2",
                                                 "Distance", "Distance-PP", "Distance-AK", "Distance-WS"])

                                for line in lines:
                                    writer.writerow(line)

                            if proc_id == master_proc:
                                print("{0} Finished writing to {1}".format(INFO, filename))

                            lines = []
                            part_counter += 1
                    touched.append((i, j))
                    touched.append((j, i))
        """
        # make sure that we flush any remaining
        # lines
        """
        if len(lines) != 0:

            if len(lines) >= line_counter:
                raise ValueError("ERROR: Detected len(lines) >={0} not written to file".format(len(lines)))

            filename = Path(outdir + "part_seqs_pairs_proc_" + str(proc_id) +"_" + str(part_counter) + ".csv")

            if filename.is_file():
                raise ValueError("ERROR: filename={0} Exists".format(filename))

            if proc_id == master_proc:

                print("{0} Writing remaining lines to {1}".format(INFO, filename))
                print("{0} Writing to {1}".format(INFO, filename))

            with open(filename, 'w', newline="\n") as fh:
                writer = csv.writer(fh, delimiter=",")
                writer.writerow(["#", "ChrSeq-1", "StartSeq-1", "EndSeq-1", "Seq-1", "HMM-State-1",
                                 "ChrSeq-2", "StartSeq-2", "EndSeq-2", "Seq-2", "HMM-State-2",
                                 "Distance", "Distance-PP", "Distance-AK", "Distance-WS"])

                for line in lines:
                    writer.writerow(line)

            if proc_id == master_proc:
                print("{0} Finished writing remaining lines to {1}".format(INFO, filename))

            #lines = []
            part_counter += 1

        error_map[proc_id] = "FINISHED SUCCESS"
    except Exception as e:
        error_map[proc_id] = str(e)


def main():

    print("{0} Starting...".format(INFO))

    configuration = read_json(filename="config_cluster.json")
    distance_type = configuration["distance_type"]
    line_counter = configuration["line_counter"]
    load = configuration["load"]

    input_file = configuration["input_file"]
    outdir = configuration["output_file"]
    num_procs = configuration["num_procs"]

    MASTER_PROC_ID = configuration["master_proc_id"]

    print("{0} Master Process is {1}".format(INFO, MASTER_PROC_ID))
    print("{0} Number of processes is {1}".format(INFO, num_procs))
    print("{0} Distance type is {1}".format(INFO, distance_type))
    print("{0} Process load {1}".format(INFO, load))

    manager = mp.Manager()
    err_map = manager.dict()

    for i in range(num_procs):
        err_map[i] = "NO_ERROR"

    procs = []

    start = configuration["start"]
    end = start + load
    for p in range(num_procs-1):
        print("{0} Process {1} works in [{2},{3})".format(INFO, p, start, end))
        procs.append(mp.Process(target=write_pair_segments_distances,
                                group=None,
                                args=(p, MASTER_PROC_ID,
                                      start, end, input_file,
                                      line_counter, outdir,
                                      distance_type, err_map)))
        procs[p].start()
        start = end
        end = end + load

    print("{0} Master Process {1} works in [{2},{3})".format(INFO, MASTER_PROC_ID, start, end))

    if "end" in configuration and configuration["end"] is not None:
        end = configuration["end"]

    # main process is working as well
    write_pair_segments_distances(proc_id=MASTER_PROC_ID,
                                  master_proc=MASTER_PROC_ID,
                                  start=start, end=end,
                                  input_file=input_file,
                                  line_counter=line_counter,
                                  outdir=outdir,
                                  distance_type=distance_type, error_map=err_map)

    # make sure we wait here
    for p in procs:
        p.join()

    for i in range(num_procs):

        if err_map[i] != "FINISHED SUCCESS":
            print("{0} ERROR detected for process {1}. Error Message {2}".format(INFO, i, err_map[i]))
        else:
            print("{0} Process {1} finished with success".format(INFO, i))
    print("{0} Finished...".format(INFO))


if __name__ == '__main__':
    main()
