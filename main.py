import pandas as pd
import matplotlib.pyplot as plt
import copy

"""
Initialize Dict based on length of dictionary and starting bit
"""
def initalize_dict(sequence):
    dictionary = {}
    # bin = string_to_binary(sequence)
    for char in sequence:
        if char not in dictionary.keys():
            dictionary[char] = len(dictionary)
    return dictionary

"""
Returns tuple(decomposition,compressed string) the compressed output 
and the decomposition of the string using LZ78 Algorithm
Note the compressed output is based on the Initialized dictionary
"""
def LZW_Compression(sequence,dictionary):
    index = 0
    decomposition = []
    compression = []
    dict = copy.deepcopy(dictionary)
    while index < len(sequence):
        temp_index = index + 1
        temp_string = sequence[index]
        while temp_index < len(sequence) and (temp_string + sequence[temp_index]) in dict.keys():
            temp_string += sequence[temp_index]
            temp_index += 1
        decomposition.append(temp_string) # Dictionary[temp_string] will give the 
        compression.append(dict[temp_string])

        if temp_index < len(sequence):
            dict[temp_string + sequence[temp_index]] = len(dict)
        index = temp_index
    return len(decomposition),compression


def string_to_binary(string_sequence):
    """Simple converter from a string sequence to a binary sequence"""
    return ''.join(format(ord(x), 'b') for x in string_sequence)

"""
Load the dataset from txt file and combine into a mega dna strand string
"""
def data_loader(filename):
     dna_dataset = filename
     df = pd.read_csv(dna_dataset, sep = "\t", names = ['sequence','class'], skiprows=1)
     combine_sequence = df['sequence'].to_list()
     join_dna = ''.join(combine_sequence)
     return join_dna

class DNASeq:

    def __init__(self,sequence,dictionary,window_size,step):
        self.sequence = sequence
        self.dictionary = initalize_dict(string_to_binary(sequence))
        self.window_size = window_size
        self.step = step
        self.lz_sequence = []
        self.normalized = []

    def split_to_window_sequences(self):
        lz = []
        start = 0
        end = self.window_size
        while end < len(self.sequence):
            trimmed_sequence = self.sequence[start:end]
            bin = string_to_binary(trimmed_sequence)
            decomp, _ = LZW_Compression(bin,self.dictionary)
            lz.append(decomp)
            start = start + self.step
            end = end + self.step
        self.lz_sequence = lz
        # Normalizes and stored as list
        self.normalize()
        return lz # Each Item is LZ78 Compressed size
    
    def normalize(self):
        self.normalized = [l/self.window_size for l in self.lz_sequence]
        return self.normalized

    def graph(self,resolution):
        plt.plot(self.normalized[0:resolution])
        plt.xlabel("x")
        plt.ylabel("Lempel-Ziv Complexity")
        plt.title("Lempel-Ziv Complexity Across DNA Sequence")



if __name__ == "__main__":
    
    join_dna = data_loader('./human.txt')
    print(f"The full length of Sequence {len(join_dna)}")

    window_size = 250
    step = 50
    dictionary = {}

    d = DNASeq(join_dna,dictionary,window_size,step)

    # Split into window sequences and apply the compression to kolmogrov Complexity Estimate
    lz = d.split_to_window_sequences()
    print(f"The number of sequence splits based on Window Size {len(lz)}")

    resolutions = [100, 1000, 10000,len(d.normalized)]
    # Get Plots from the normalized sequence for different resolution measures
    for res in resolutions:
        d.graph(res)
        # plt.show()
        plt.savefig(f'{res}.png', dpi = 600)
    