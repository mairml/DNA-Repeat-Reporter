import sys, re #used regular expression findinter

header = [] #input FASTA file header
motifs = {} #dictionary of all repeated sequences
indices = {} #index of each base in a repeated sequence
i = 0

#Parse an input FASTA file
def read_file(input_file):
    print("Reading sequence file...")
    file1 = open(input_file, 'r')
    seq = [] #list for sequence lines, helps for reading a multi-line FASTA file
    for line in file1:
        # parsing the file
        line = line.rstrip('\n') # getting rid of new line characters
        if line.startswith('>'):
            header.append(line.strip('\n').strip('>'))
        else:
            line = line.strip('\n')
            seq.append(line)
    sequence = "".join(seq) #take all the seq lines in a multi-line FASTA file, put them all in the same string
    seq = []#empty the list because we don't need it any longer
    return sequence

#Take index of each nucleotide in a repeat and stores them in dictionary, helps generate alignment file.
def index(window, repeat_indices):
    for i in range(0,len(repeat_indices)):
        for x in range(0,len(window)):
            base = window[x]
            index = repeat_indices[i]+x
            indices[index]=base
    return(indices)

#Finds the indexrange of a repeated sequence and converts format for printing
def repeat_range(repeat_indices, dynamic_window):
    expect = None
    run = []
    result = [run]
    #Go through indices in list and group them if they are consecutive
    for index in repeat_indices:
        if (expect is None) or (index == expect):
            run.append(index)
        else:
            run = [index]
            result.append(run)
        expect = index + 1
    final = []
    #Format a range for printing
    for items in result:
        if len(items)>1:
            #NOTE: Need to add 1 because Python counts from 0
            first = str(items.pop(0)+1)
            last = str(items.pop(-1)+dynamic_window)
            final.append([first + "-" + last])
        elif len(items)==1:
            items = list(map(lambda x: x + 1, items))
            first = str(items.pop(0))
            last = str(int(first)-1+dynamic_window)
            final.append([first + "-" + last])
        else:
            final.append("ERROR") #Place holder that should not be used...make sure data formatting worked.
    printout = str(final).replace('[', '').replace(']','').replace("'",'') #turn list into easy-to-read text
    return printout

#Find repeated sequence and its index
def repeated(window, i,dynamic_window):
    if sequence.count(window) > 1:
        size = len(window)
        repeat_indices = [base.start() for base in re.finditer(window, sequence)] #generate list of indices where each sequence occurs'
        printout = repeat_range(repeat_indices,dynamic_window)#repeat_range function defined above
        repeat_index = []
        index(window, repeat_indices)#index function defined above
        repeat_indices= list(map(lambda x: x+1, repeat_indices))
        for num in repeat_indices:
            new = str(num) + "-" + str((num-1)+size)
            repeat_index.append(new)
        if window not in motifs:
            motifs[window] = repeat_indices
            report_outfile.write("####################################################\n" +
                                 "#REPEAT of length " + str(size) + " found!" + "\n" +
                                 "#Repeated sequence found at index: " + printout + "\n" +
                                 "#Repeat Identity" + "\n" +
                                 "#" + window + "\n")
        return(window)

#Find and take the index of reverse complement sequences
def revComplement(window, i, window_size, dynamic_window):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    revComp = "".join(complements.get(base, base) for base in reversed(window))  # create reverse complement list, turn it into a string
    if revComp in sequence:
        repeat_index = sequence.find(revComp)
        if repeat_index < i + dynamic_window and i < repeat_index + dynamic_window:  # ensure the repeated sequence does not contain a fragment of the original window (no mathematical overlap of index ranges)
            pass
        else:
            revComp_re = revComp
            repeat_indices = [base.start() for base in re.finditer(revComp_re, sequence)]#Take the index for output sequence alignment
            printout = repeat_range(repeat_indices, dynamic_window) #Format for printing
            index(revComp, repeat_indices)
            if sequence.count(window) > 1:
                if revComp not in motifs:
                    motifs[window] = repeat_indices
                    report_outfile.write("#Reverse complement sequence found at index: " + printout + "\n"+
                                         "#Reverse Complement Sequence Identity" + "\n" +
                                         "#" + revComp + "\n")
                    return "continue" #repeat found, keep searching for repeats in dynamic window
            else:
                if revComp not in motifs:
                    motifs[window] = repeat_indices
                    index(window, [i])
                    length = str(len(window))
                    printout2 = repeat_range([i],dynamic_window)
                    report_outfile.write("####################################################\n" + "#REPEAT of length " + length + " found! \n" +
                                         "#Repeated sequence found at index: " + printout2 + "\n" +
                                         "#Sequence IDENTITY" + "\n" +"#" + window + "\n"+
                                         "#Reverse complement sequence found at index: " + printout + "\n" +
                                         "#Reverse Complement Sequence Identity" + "\n" +
                                         "#" + revComp + "\n")
                    return "continue" #repeat found, keep searching for repeats in dynamic window
    else:
        if sequence.count(window) <=1:
            return "stop" #return that there's no longer a repeat, so we should stop


#Main code... parse through sequence with sliding window

print("Running the Repeat Finder")
#user-defined input file
infile = sys.argv[1]

#user-defined window-size
input_window = int(sys.argv[2])

#user-defined report file name
report_outfile = open(sys.argv[3],'w')

#user-defined alignment file name
seq_outfile = open(sys.argv[4], 'w')
sequence = read_file(infile)


print("Finding repeat sequences and creating report. This may take a while...")
window_size = input_window
report_outfile.write("Repeats of minimum length " + str(window_size) + " in " + str(header).strip('[').strip(']').strip("'") + "\n")


def loop(sequence,window_size,):
    for i in range(0,len(sequence),+1): #generate sliding window
        window = sequence[i:i+window_size]
        #The following is an optional progress statement. Finding repeats takes some time, so this prints out alerts you the program is still running
        print("Parsing nucleotide " + str(i) +"/" + str(len(sequence)) + "\r")
        #Find a repeat of  the window:
        if window not in motifs: #make sure the repeat has not already been identified and documented (prevents indentical windows)
            if len(window) == window_size: #ensure the window is the right size
                #find repeats that are multiples of the window size
                window = sequence[i:i + window_size]
                repeated(window, i,window_size)
                revComplement(window, i, window_size,window_size)
                #find repeats that are not multiples of the window size (i.e., window_size =3, find repeats that are 4 and 5 bp long)
                for z in range(1,window_size):
                    if (window_size+z)%window_size != 0:
                        dynamic_window = window_size+z #increase window size
                        window = sequence[i:i + window_size+z]
                        repeated(window,i,dynamic_window)
                        n = revComplement(window, i, window_size,dynamic_window)
                        if n == "stop": #no point in continuing the count if repeat of given size isn't found (i.e., window_size = 3 (AAA), no point in seeing if repeats of length 5 (AAATT) exist if length 4 (AAAT) doesn't exist)
                            break
                            
loop(sequence,window_size)

#Create output FASTA file to visualize repeats through sequence alignment
print("Report complete. Generating file for sequence alignment...")
new_sequence=[] #Build list of repeated sequence bases or unknown nucleotide N where no repeated sequence occurs
for i in range(0,len(sequence)):
    if i in indices:
        new_sequence.append(indices[i])
    else:
        new_sequence.append("N")

output_sequence = "".join(new_sequence)#put the list into a string
#Print output FASTA file
seq_outfile.write(">REPEATS "+str(header).strip('[').strip(']').strip("'") +"\n")
seq_outfile.write(output_sequence)
print("Program successfully completed.")
