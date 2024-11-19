import sys

##Step 1: Argument Parsing and File Selection.
if len(sys.argv) == 3:
    fasta_file = sys.argv[1]
    enzyme_file = sys.argv[2]
else:
    #Prompting the user to enter file paths if not provided as arguments.
    fasta_file = input("Please enter the path of the FASTA file containing the nucleotide sequence: ")
    enzyme_file = input("Please enter the path of the file containing the restriction enzymes: ")

##Step 2: Initializing Variables and Reading the Sequence Data.
sequence_name = "Unknown"
nucleotide_seq = ""

#Reading the sequence file.
seq_file = open(fasta_file, 'r')
first_line = seq_file.readline().strip()

#Assuming the first line is always a header (starts with ">").
sequence_name = first_line[1:]  #Removing the ">" and taking the rest as the name.

#Continuing to read the rest of the sequence lines (after the header).
for line in seq_file:
    nucleotide_seq += line.strip()  #Adding each line of the sequence.
seq_file.close()  #Closing the file.

##Step 3: Calculating Sequence Length and Outputting Initial Information.
sequence_length = len(nucleotide_seq)

#Printing the formatted header information.
print("Restriction enzyme analysis of sequence from file " + fasta_file + ".")
print("Cutting with enzymes found in file " + enzyme_file + ".")
print("---------------------------------------------------------------")
print("Sequence name: " + sequence_name)
print("Sequence is " + str(sequence_length) + " bases long.")
print("---------------------------------------------------------------")

##Step 4: Analyzing Each Enzyme for Cutting Sites in the Sequence.
enzyme_file = open(enzyme_file, 'r')
for line in enzyme_file:
    #Extracting enzyme name and recognition sequence.
    enzyme_data = line.strip().split(';')
    enzyme_name = enzyme_data[0]
    recognition_seq = enzyme_data[1]  #Assuming the cleavage marker is always "^".

    #Locating the cleavage position in the recognition sequence.
    cleavage_pos = recognition_seq.index("^")
    recognition_seq = recognition_seq.replace("^", "")

    #Finding all cutting sites.
    cutting_sites = []
    start = nucleotide_seq.find(recognition_seq)  #Initial search for the recognition sequence.

    #Looping to find each occurrence of the recognition sequence.
    while start != -1:
        cutting_sites.append(start + cleavage_pos)  #Position of the actual cut.
        start = nucleotide_seq.find(recognition_seq, start + 1)  #Searching for the next occurrence.

    ##Step 5: Generating and Outputting Fragments Based on Cutting Sites.
    if cutting_sites:
        #Including the start (0) and end (sequence length) positions in the cutting sites.
        all_positions = [0] + cutting_sites + [len(nucleotide_seq)]
        fragments = []  #Initializing an empty list to store fragments.

        #Looping through positions in pairs to get each fragment.
        for i in range(len(all_positions) - 1):
            start_pos = all_positions[i]
            end_pos = all_positions[i + 1]
            fragment = nucleotide_seq[start_pos:end_pos]  #Extracting the fragment from the sequence.
            fragments.append(fragment)  #Adding the fragment to the fragments list.

        #Printing information about cutting sites and fragments.
        print("There are " + str(len(cutting_sites)) + " cutting sites for " + enzyme_name +
              ", cutting at " + recognition_seq[:cleavage_pos] + "^" + recognition_seq[cleavage_pos:] + ".")
        print("There are " + str(len(fragments)) + " fragments:\n")

        #Simplified formatting for each fragment output.
        position = 1  #Starting position for each fragment.
        for fragment in fragments:
            fragment_length = len(fragment)
            print("Length- " + str(fragment_length))

            #Printing each fragment in lines of 60 bases, grouped in sets of 10.
            index = 0
            while index < fragment_length:
                line = fragment[index:index+60]  #Getting the next 60 bases.
                formatted_line = ""

                #Looping to group each line into sets of 10 bases.
                for group_start in range(0, len(line), 10):
                    formatted_line += line[group_start:group_start+10] + " "

                #Removing the extra space at the end and printing with position.
                formatted_line = formatted_line.strip()
                print(str(position) + "\t" + formatted_line)

                #Updating position and index.
                position += len(line)
                index += 60

            print()  #Blank line between fragments.

        print("---------------------------------------------------------------")
    else:
        print("There are no sites for " + enzyme_name + ".")
        print("---------------------------------------------------------------")
enzyme_file.close()  #Closing the enzyme file.
