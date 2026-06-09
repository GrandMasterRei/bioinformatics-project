from ahmet_bioproject.muscle_core import my_muscle_msa

test_seqs = ["ATGC", "ATGCGT", "AGCGT"]
for seq in my_muscle_msa(test_seqs):
    print(seq)