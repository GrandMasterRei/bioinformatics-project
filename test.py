
from ahmet_bioproject.muscle_core import my_muscle_msa

# Test etmek istediğimiz örnek DNA/RNA dizileri
ornek_diziler = ["ATGC", "ATGCGT", "AGCGT"]

print("Algoritma çalıştırılıyor...")
sonuclar = my_muscle_msa(ornek_diziler)

print("\nHizalanmış Diziler (Çıktı):")
for s in sonuclar:
    print(s)