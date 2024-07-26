# This is our UBE2C protein sequence
UBE2C_sequence = "MASQNRDPAATSVAAARKGAEPSGGAARGPVGKRLQQELMTLMMSGDKGISAFPESDNLF" \
                  "KWVGTIHGAAGTVYEDLRYKLSLEFPSGYPYNAPTVKFLTPCYHPNVDTQGNICLDILKE" \
                  "KWSALYDVRTILLSIQSLLGEPNIDSPLNTHAAELWKNPTAFKKYLQETYSKQVTSQEP"
# SSE values
true_structure = [
    ("Helix", 30, 45),
    ("Beta strand", 50, 54),
    ("Beta strand", 61, 68),
    ("Beta strand", 78, 84),
    ("Turn", 87, 91),
    ("Beta strand", 95, 100),
    ("Beta strand", 111, 113),
    ("Helix", 116, 118),
    ("Turn", 119, 121),
    ("Helix", 128, 140),
    ("Helix", 150, 155),
    ("Helix", 159, 172)
]
true_labels = ["U"] * len(UBE2C_sequence)
# I'm renaming the structures' names at SSE again for comparison.
for label, start, end in true_structure:
    if label == "Beta strand":
        true_labels[start - 1:end] = ["E"] * (end - start + 1)
    else:
        true_labels[start - 1:end] = [label[0]] * (end - start + 1)
# Our parameters to do in Chou-Fasman method.
window_size = 6
helix_threshold = 1
sheet_threshold = 1
turn_threshold = 0.000075

# Propensity table. I've accepted some letters as if they would represent those that begin
# with the same initial letter. ( Like Glycine's values for G )
propensity_table = {
    "A": [1.42, 0.83, 0.66, 0.06, 0.076, 0.035, 0.058],
    "R": [0.98, 0.93, 0.95, 0.070, 0.106, 0.099, 0.085],
    "N": [0.67, 0.89, 1.56, 0.161, 0.083, 0.191, 0.091],
    "D": [1.01, 0.54, 1.46, 0.147, 0.110, 0.179, 0.081],
    "C": [0.70, 1.19, 1.19, 0.149, 0.050, 0.117, 0.128],
    "Q": [1.11, 1.10, 0.98, 0.074, 0.098, 0.037, 0.098],
    "E": [1.39, 1.17, 0.74, 0.056, 0.060, 0.077, 0.064],
    "G": [0.57, 0.75, 1.56, 0.102, 0.085, 0.190, 0.152],
    "H": [1.00, 0.87, 0.95, 0.140, 0.047, 0.093, 0.054],
    "I": [1.08, 1.60, 0.47, 0.043, 0.034, 0.013, 0.056],
    "L": [1.41, 1.30, 0.59, 0.061, 0.025, 0.036, 0.070],
    "K": [1.14, 0.74, 1.01, 0.055, 0.115, 0.072, 0.095],
    "M": [1.45, 1.05, 0.60, 0.068, 0.082, 0.014, 0.055],
    "F": [1.13, 1.38, 0.60, 0.059, 0.041, 0.065, 0.065],
    "P": [0.57, 0.55, 1.52, 0.102, 0.301, 0.034, 0.068],
    "S": [0.77, 0.75, 1.43, 0.120, 0.139, 0.125, 0.106],
    "T": [0.83, 1.19, 0.96, 0.086, 0.108, 0.065, 0.079],
    "W": [1.08, 1.37, 0.96, 0.077, 0.013, 0.064, 0.167],
    "Y": [0.69, 1.47, 1.14, 0.082, 0.065, 0.114, 0.125],
    "V": [1.06, 1.70, 0.50, 0.062, 0.048, 0.028, 0.053]
}

def is_helix(subsequence):
    return sum(propensity_table[aa][0] > helix_threshold for aa in subsequence) >= 4

def is_sheet(subsequence):
    return sum(propensity_table[aa][1] > sheet_threshold for aa in subsequence) >= 4

def is_turn(subsequence):
    bend_score = propensity_table[subsequence[2]][3] * propensity_table[subsequence[3]][4] * \
                 propensity_table[subsequence[4]][5] * propensity_table[subsequence[5]][6]
    return bend_score >     turn_threshold

# Looking for overlaps here.
def check_overlap(label_seq, start, end):
    for i in range(max(0, start), min(len(label_seq), end)):
        if label_seq[i] != 'U':
            return True
    return False

# Chou-Fasman method processing here.
predicted_labels = ["U"] * len(UBE2C_sequence)

for i in range(len(UBE2C_sequence) - window_size + 1):
    subsequence = UBE2C_sequence[i:i + window_size]

    # Alpha-helix checking here.
    if is_helix(subsequence):
        # Extending.
        for j in range(i, i + window_size):
            if not check_overlap(predicted_labels, j, j + window_size):
                predicted_labels[j] = "H"

    # Same processes for beta strand.
    if is_sheet(subsequence):
        for j in range(i, i + window_size):
            if not check_overlap(predicted_labels, j, j + window_size):
                predicted_labels[j] = "E"

    # Calculating both of them a and b for checking overlap and decide.
    if is_helix(subsequence) and is_sheet(subsequence):
        overlap_start = max(0, i)
        overlap_end = min(len(predicted_labels), i + window_size)
        sp_a = sum(propensity_table[aa][0] for aa in predicted_labels[overlap_start:overlap_end])
        sp_b = sum(propensity_table[aa][1] for aa in predicted_labels[overlap_start:overlap_end])

        if sp_a > sp_b:
            for j in range(overlap_start, overlap_end):
                predicted_labels[j] = "H"
        else:
            for j in range(overlap_start, overlap_end):
                predicted_labels[j] = "E"

    # Turn checking here.
    if is_turn(subsequence):
        for j in range(i, i + 4):
            if not check_overlap(predicted_labels, j, j + 4):
                predicted_labels[j] = "T"


# Confusion matrix.
confusion_matrix = {
    "H": {"H": 0, "E": 0, "T": 0, "U": 0},
    "E": {"H": 0, "E": 0, "T": 0, "U": 0},
    "T": {"H": 0, "E": 0, "T": 0, "U": 0},
    "U": {"H": 0, "E": 0, "T": 0, "U": 0},
}

for true_label, start, end in true_structure:
    true_label_seq = true_labels[start - 1:end]
    predicted_label_seq = predicted_labels[start - 1:end]

    for true, pred in zip(true_label_seq, predicted_label_seq):
        confusion_matrix[true][pred] += 1

# Metrics.
def calculate_metrics(matrix, label):
    tp = matrix[label][label]
    fp = sum(matrix[pred_label][label] for pred_label in matrix if pred_label != label)
    fn = sum(matrix[label][true_label] for true_label in matrix if true_label != label)
    tn = sum(matrix[pred_label][true_label] for pred_label in matrix
             for true_label in matrix if true_label != label and pred_label != label)

    total_elements = sum(value for inner_dict in matrix.values() for value in inner_dict.values())

    precision = tp / (tp + fp) if tp + fp != 0 else 0
    recall = tp / (tp + fn) if tp + fn != 0 else 0
    accuracy = (tp + tn) / total_elements if total_elements != 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall != 0 else 0

    return precision, recall, accuracy, f1_score



precision_H, recall_H, accuracy_H, f1_score_H = calculate_metrics(confusion_matrix, "H")
precision_E, recall_E, accuracy_E, f1_score_E = calculate_metrics(confusion_matrix, "E")
precision_T, recall_T, accuracy_T, f1_score_T = calculate_metrics(confusion_matrix, "T")
precision_U, recall_U, accuracy_U, f1_score_U = calculate_metrics(confusion_matrix, "U")

# Printing the confusion matrix and scores.
print("Confusion Matrix:")
print("Predicted\tTrue\tH\tE\tT\tU")
for true_label in confusion_matrix:
    print(f"{true_label}\t\t", end="")
    for pred_label in confusion_matrix[true_label]:
        print(f"{confusion_matrix[true_label][pred_label]}\t", end="")
    print()

print("\nMetrics:")
print("Label\t\tPrecision\tRecall\t\tAccuracy\tF1-Score")
print(f"H\t\t{precision_H:.4f}\t\t{recall_H:.4f}\t\t{accuracy_H:.4f}\t\t{f1_score_H:.4f}")
print(f"E\t\t{precision_E:.4f}\t\t{recall_E:.4f}\t\t{accuracy_E:.4f}\t\t{f1_score_E:.4f}")
print(f"T\t\t{precision_T:.4f}\t\t{recall_T:.4f}\t\t{accuracy_T:.4f}\t\t{f1_score_T:.4f}")
print(f"U\t\t{precision_U:.4f}\t\t{recall_U:.4f}\t\t{accuracy_U:.4f}\t\t{f1_score_U:.4f}")
