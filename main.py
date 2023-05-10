import pandas as pd
import matplotlib.pyplot as plt


"""
Initialize Dict based on length of dictionary and starting bit
"""
def initalize_dict(sequence):
    Dictionary = {}
    for char in sequence:
        if char not in Dictionary.keys():
            Dictionary[char] = len(Dictionary)
    return Dictionary

"""
Returns tuple(decomposition,compressed string) the compressed output 
and the decomposition of the string
Note the compressed output is based on the Initialized dictionary
"""
def LZW_Compression(sequence):
    index = 0
    decomposition = []
    compression = []
    Dictionary = initalize_dict(sequence)
    # print(Dictionary)
    while index < len(sequence):

        temp_index = index + 1
        temp_string = sequence[index]
        # substring = temp_string + sequence[temp_index]

        while temp_index < len(sequence) and (temp_string + sequence[temp_index]) in Dictionary.keys():
            temp_string += sequence[temp_index]
            temp_index += 1
            # try:
            #     substring = temp_string + sequence[temp_index]
            # except:
            #     continue

        decomposition.append(temp_string) # Dictionary[temp_string] will give the 
        compression.append(Dictionary[temp_string])

        if temp_index < len(sequence):
            Dictionary[temp_string + sequence[temp_index]] = len(Dictionary)

        index = temp_index

    return len(decomposition),compression

def string_to_binary(string_sequence):
    """Simple converter from a string sequence to a binary sequence"""
    return ''.join(format(ord(x), 'b') for x in string_sequence)

# print(string_to_binary(''))

"""Generator that yields the window of the sequence"""
def window(sequence,window_size = 250,step = 50):
    start = 0
    end = window_size
    while end < len(sequence):
        yield sequence[start:end]
        start = start + step
        end = end + step

def lz_sequence(sequence):
    lz = []
    for w in window(sequence):
        bin = string_to_binary(w)
        decomp, _ = LZW_Compression(bin)
        lz.append(decomp)
    return lz

def normalize_lz(lz_sequence,window_size):
    return [l/window_size for l in lz_sequence]

def graph(lz_sequence):
    plt.plot(lz_sequence)
    plt.xlabel("Position")
    plt.ylabel("Lempel-Ziv Complexity")
    plt.title("Lempel-Ziv Complexity Across DNA Sequence")

if __name__ == "__main__":
    dna_dataset = './human.txt'
    df = pd.read_csv(dna_dataset, sep = "\t", names = ['sequence','class'], skiprows=1)
    combine_sequence = df['sequence'].to_list()[:5]
    join_dna = ''.join(combine_sequence)
    print(len(join_dna))
    lz  = lz_sequence(join_dna)
    print(lz)
   
    n_lz_sequence = normalize_lz(lz, 250)
    print(n_lz_sequence)
    # resolutions = [100, 1000, 10000]

    # for res in resolutions:
    #     graph(n_lz_sequence[0:res])
    #     plt.savefig(f'{res}.png', dpi = 600)
    
    # graph(n_lz_sequence)
    # plt.savefig(f'full.png', dpi = 600)

    # print(df.loc[:,'sequence'])
    # print(df)
    # s = '1001111011000010'
    # Sample sequence of LZW
    # 1 0 10
    # 0 0 00
    # 0 1 01
    # 1 1 11
    # 11 1 111
    # 1 0 0 1 11 10 11 00 00 10
    # 0 1 1 0 5 2 5 3 3 2
    # print(LZW_Compression(s))



