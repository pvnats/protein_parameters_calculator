from tkinter import *
from tkinter import filedialog

# Protein parameters calculation by pvnats
# Based on ProtParam https://www.expasy.org/resources/protparam
# Computes amino acid composition, atomic composition, molecular weight,
# extinction coefficients and theoretical pI from a protein sequence

# Amino acids dictionary
# One letter : [3 letters, Mw, C, H, N, O, S, pKa]
amino_acids_dict = {
    "A": ["Ala", 89.09404, 3, 5, 1, 1, 0],
    "C": ["Cys", 121.15404, 3, 5, 1, 1, 1, 9],
    "D": ["Asp", 133.10384, 4, 5, 1, 3, 0, 4.05],
    "E": ["Glu", 147.13074, 5, 7, 1, 3, 0, 4.45],
    "F": ["Phe", 165.19184, 9, 9, 1, 1, 0],
    "G": ["Gly", 75.06714, 2, 3, 1, 1, 0],
    "H": ["His", 155.15634, 6, 7, 3, 1, 0, 5.98],
    "I": ["Ile", 131.17464, 6, 11, 1, 1, 0],
    "K": ["Lys", 146.18934, 6, 12, 2, 1, 0, 10],
    "L": ["Leu", 131.17464, 6, 11, 1, 1, 0],
    "M": ["Met", 149.20784, 5, 9, 1, 1, 1],
    "N": ["Asn", 132.11904, 4, 6, 2, 2, 0],
    "P": ["Pro", 115.13194, 5, 7, 1, 1, 0],
    "Q": ["Gln", 146.14594, 5, 8, 2, 2, 0],
    "R": ["Arg", 174.20274, 6, 12, 4, 1, 0, 12],
    "S": ["Ser", 105.09344, 3, 5, 1, 2, 0],
    "T": ["Thr", 119.12034, 4, 7, 1, 2, 0],
    "V": ["Val", 117.14784, 5, 9, 1, 1, 0],
    "W": ["Trp", 204.22844, 11, 10, 2, 1, 0],
    "Y": ["Tyr", 181.19124, 9, 9, 1, 2, 0, 10],
}


def make_sequence(my_sequence):  # Making the sequence list of the analyzing protein
    sequence_list = []
    for s in my_sequence:
        if s.upper() in amino_acids_dict.keys():
            sequence_list.append(s.upper())
            
    result_text.insert("end", f"Your sequence:\n")
    
    # Sequence chunking        
    chunked_list = [sequence_list[i:i + 10] for i in range(0, len(sequence_list), 10)]
    counter = 0
    for i in range(0, len(chunked_list), 5):  # Chunked sequence printing
        for chunk in chunked_list[i:i + 5]:
            counter += len(chunk)
            result_text.insert("end", f"{''.join(chunk)} ")
        result_text.insert("end", f' {counter}\n')
    return sequence_list

def count_aa(sequence_list):  # Amino acids counting
    total_aa_number = len(sequence_list)

    temp_aa_dict = {}
    for aa in sequence_list:
        temp_aa_dict[aa] = temp_aa_dict.get(aa, 0) + 1

    result_text.insert("end", f"\nAmino acid composition: \n")

    for key, value in sorted(temp_aa_dict.items()):
        result_text.insert(
            "end",
            f"{amino_acids_dict[key][0]} ({key}) = {value}\t"
            f"\t{100 * value / total_aa_number:.1f}% \n",
        )

    result_text.insert("end", f"\nNumber of amino acids: {total_aa_number}\n")
    return temp_aa_dict


def count_aa_atoms(my_sequence_dict):  # Amino acid composition calculation
    atom_c, atom_h, atom_n, atom_o, atom_s = 0, 0, 0, 0, 0

    for key, value in my_sequence_dict.items():
        atom_c += amino_acids_dict[key][2] * value
        atom_h += amino_acids_dict[key][3] * value
        atom_n += amino_acids_dict[key][4] * value
        atom_o += amino_acids_dict[key][5] * value
        atom_s += amino_acids_dict[key][6] * value

    result_text.insert("end", f"\nAtomic composition:")
    result_text.insert(
        "end",
        f"\nCarbon  \tC\t{atom_c}"
        f"\nHydrogen\tH\t{atom_h + 2}"
        f"\nNitrogen\tN\t{atom_n}"
        f"\nOxygen  \tO\t{atom_o + 1}"
        f"\nSulfur  \tS\t{atom_s}\n",
    )

    result_text.insert(
        "end", f"\nFormula: C{atom_c}H{atom_h + 2}N{atom_n}O{atom_o + 1}S{atom_s}\n"
    )
    result_text.insert(
        "end",
        f"Total number of atoms: {atom_c + atom_h + atom_n + atom_o + atom_s + 3}\n",
    )


def calculate_mol_weight(my_sequence_dict):  # Molecular weight calculation
    water_weight = 18.0153
    mol_weight = 0

    for key, value in my_sequence_dict.items():
        mol_weight += (amino_acids_dict[key][1] - water_weight) * value
    result_text.insert("end", f"\nMolecular weight: {mol_weight + water_weight:.2f}\n")
    return mol_weight


def calculate_ext_coef(
    my_sequence_dict, mol_weight
):  # Extinction coefficients calculation
    ext_coef_cformed = (
        my_sequence_dict.get("Y", 0) * 1490
        + my_sequence_dict.get("W", 0) * 5500
        + my_sequence_dict.get("C", 0) // 2 * 125
    )
    ext_coef_reduced = (
        my_sequence_dict.get("Y", 0) * 1490 + my_sequence_dict.get("W", 0) * 5500
    )
    abs_cformed = ext_coef_cformed / mol_weight
    abs_reduced = ext_coef_reduced / mol_weight

    result_text.insert("end", f"\nExt. coefficient (M-1 cm-1): {ext_coef_cformed}\n")
    result_text.insert(
        "end",
        f"Abs (0.1%): {abs_cformed:.3f}, assuming all pairs of Cys residues form cystines\n",
    )

    result_text.insert("end", f"\nExt. coefficient (M-1 cm-1): {ext_coef_reduced}\n")
    result_text.insert(
        "end", f"Abs (0.1%): {abs_reduced:.3f}, assuming all Cys residues are reduced\n"
    )


def calculate_pI(my_sequence_dict, sequence_list):
    """Protein isoelectric point definition using Henderson-Hasselbach equation
    http://isoelectric.org/www_old/files/isoelectric-point-theory.html

    pKa dictionaries according Protparam:
    Bjellqvist, B., Basse, B., Olsen, E., & Celis, J. E. (1994).
    Reference points for comparisons of two-dimensional maps of proteins from
    different human cell types defined in a pH scale where isoelectric points
    correlate with polypeptide compositions. Electrophoresis, 15(1), 529–539.
    doi:10.1002/elps.1150150171
    Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
    Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F.
    The focusing positions of polypeptides in immobilized pH gradients can be
    predicted from their amino acid sequences. Electrophoresis 1993, 14,
    1023-1031.

    Calculations based on the examplary program proposed on
    http://isoelectric.org/www_old/files/practise-isoelectric-point.html"""

    nter_dict = {
        "Nter": 7.5,
        "A": 7.59,
        "M": 7.00,
        "S": 6.93,
        "P": 8.36,
        "T": 6.82,
        "V": 7.44,
        "E": 7.70,
    }  # N-ter pKa
    cter_dict = {"Cter": 3.55, "D": 4.55, "E": 4.75}  # C-Ter pKa

    # N-terminal pKa of N-ter AA of my_sequence
    if sequence_list[0] in nter_dict.keys():
        pKa_nter = nter_dict[sequence_list[0]]
    else:
        pKa_nter = nter_dict["Nter"]

    # C-terminal pKa of C-ter AA of my_sequence
    if sequence_list[-1] in cter_dict.keys():
        pKa_cter = cter_dict[sequence_list[-1]]
    else:
        pKa_cter = cter_dict["Cter"]

    ph = 6.50  # starting point pI = 6.5
    net_charge = 0.0  # net charge in given pH

    ph_prev = 0.0  # 0-14 is possible pH range
    ph_next = 14.0
    temp_ph = 0.0
    EPSILON = 0.01  # defined precision

    # Calculations using bisection
    while True:
        ch_cter = -1 / (1 + pow(10, (pKa_cter - ph)))
        ch_d = -my_sequence_dict.get("D", 0) / (
            1 + pow(10, (amino_acids_dict["D"][7] - ph))
        )
        ch_e = -my_sequence_dict.get("E", 0) / (
            1 + pow(10, (amino_acids_dict["E"][7] - ph))
        )
        ch_c = -my_sequence_dict.get("C", 0) / (
            1 + pow(10, (amino_acids_dict["C"][7] - ph))
        )
        ch_y = -my_sequence_dict.get("Y", 0) / (
            1 + pow(10, (amino_acids_dict["Y"][7] - ph))
        )
        ch_nter = 1 / (1 + pow(10, (ph - pKa_nter)))
        ch_h = my_sequence_dict.get("H", 0) / (
            1 + pow(10, (ph - amino_acids_dict["H"][7]))
        )
        ch_k = my_sequence_dict.get("K", 0) / (
            1 + pow(10, (ph - amino_acids_dict["K"][7]))
        )
        ch_r = my_sequence_dict.get("R", 0) / (
            1 + pow(10, (ph - amino_acids_dict["R"][7]))
        )
        net_charge = ch_cter + ch_d + ch_e + ch_c + ch_y + ch_nter + ch_h + ch_k + ch_r
        temp_ph = ph

        if net_charge < -EPSILON:  # out of range, pH value must be smaller
            ph = ph - ((ph - ph_prev) / 2)
            ph_next = temp_ph
        elif net_charge > 0:  # to small pH value, have to be increased
            ph = ph + ((ph_next - ph) / 2)
            ph_prev = temp_ph
        else:
            break

    result_text.insert("end", f"\nTheoretical pI: {ph:.2f}\n\n")


def start_calculation():  # Start    
    result_text.delete("1.0", "end")  # Cleaning the result window
    result_window.deiconify()  
    my_sequence = prot_seq_text.get("1.0", "end")
    sequence_list = make_sequence(my_sequence)
    my_sequence_dict = count_aa(sequence_list)
    count_aa_atoms(my_sequence_dict)
    mol_weight = calculate_mol_weight(my_sequence_dict)
    calculate_ext_coef(my_sequence_dict, mol_weight)
    calculate_pI(my_sequence_dict, sequence_list)
    btn_save.config(state="normal", cursor="hand2")  # Activating the Save as.. button
    input_on_left_click_back()

def input_on_left_click(event):  # Preparing the sequence window for the new data
    global prot_seq_text_click_flag
    if prot_seq_text_click_flag == False:
        prot_seq_text.delete("1.0", "end")  # Cleaning the sequence window 
        prot_seq_text.config(fg="black")  # Changing the font color
        prot_seq_text_click_flag = True  # Stop the window mouse click event
        
def input_on_left_click_back():  # to back the initial prot_seq_text window settings
    global prot_seq_text_click_flag
    prot_seq_text_click_flag = False  # Start the window mouse click event
    prot_seq_text.config(fg="grey")
    
def save_file():
    filepath = filedialog.asksaveasfilename(
        defaultextension="txt",
        initialfile="new_seq.txt",
        filetypes=[("Text files", ".txt")],
    )
    if filepath != "":
        text = result_text.get("1.0", "end")
        with open(filepath, "w") as file:
            file.write(text)
            
def on_closing_result_window():  # Result window could not be closed
    pass

def copy_paste(e):  # Fixing the hotkeys copy-paste problem when you use RU keyboard
    if e.keycode == 86 and e.keysym != 'v':
        e.widget.event_generate('<<Paste>>')
    elif e.keycode == 67 and e.keysym != 'c':
        e.widget.event_generate('<<Copy>>')
    elif e.keycode == 88 and e.keysym != 'x':
        e.widget.event_generate('<<Cut>>')
        
        
root = Tk()
root.geometry("600x400+800+300")  # Main window
root.resizable(False, False)
root.title("Tk_Protein_Parameters")

frame_title = Frame(root, width=400, height=50)
frame_title.pack()

title_l = Label(
    frame_title,
    text="Protein parameters calculator",
    font=("Arial", 16),
    height=2,
)
title_l.pack()

frame_top = Frame(root, width=600, height=50)
frame_top.pack()

# 'Calculate' button
btn_calc = Button(
    frame_top,
    text="Calculate",
    font=("Arial", 12),
    bg="lightgrey",
    cursor="hand2",
    command=start_calculation,
)
btn_calc.place(width=150, height=30, relx=0.36, rely=0.2)

# 'Save as..' button
btn_save = Button(
    frame_top,
    text="Save result",
    font=("Arial", 12),
    bg="lightgrey",
    state="disabled",
    command=save_file,
)
btn_save.place(width=100, height=30, relx=0.76, rely=0.2)

frame_bottom = Frame(root, width=600, height=300)
frame_bottom.pack()

# Sequence input 
prot_seq_text = Text(frame_bottom, wrap="char", width=72, height=17, fg="grey")
prot_seq_text.place(x=10, y=10)
prot_seq_text.insert(
    "1.0", " Your protein sequence: (Example: SKYAVK...)"
)
# prot_seq_text_click_flag become 'False' before first seq input and after every job
# to control input_on_left_click event and let to edit your seq before calculations
prot_seq_text_click_flag = False  
prot_seq_text.bind("<Button-1>", input_on_left_click)
prot_seq_text.bind("<Control-Key>", copy_paste)

# Result window
result_window = Toplevel() 
result_window.title("Results")
result_window.geometry("530x570+350+250")
result_window.protocol("WM_DELETE_WINDOW", on_closing_result_window)
result_window.withdraw()  # Hide the result window

result_text = Text(result_window, wrap="word", width=45, highlightthickness=6, padx=10, pady=5)
result_text.pack(
    anchor=N + W, 
    fill=BOTH, 
    side=LEFT, 
    expand=True
    )
scroll_text = Scrollbar(
    result_window, width=20, highlightbackground="green", command=result_text.yview
)
scroll_text.pack(fill=Y, side=LEFT, expand=False)
result_text.config(yscrollcommand=scroll_text.set)


root.mainloop()