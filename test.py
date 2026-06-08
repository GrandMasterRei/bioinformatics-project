from ahmet_bioproject import my_muscle_msa

test_seqs = ["ATGC", "ATGCGT", "AGCGT"]
sonuclar = my_muscle_msa(test_seqs)

for seq in sonuclar:
    print(seq)