def get_kmers(seq, k=2):

    kmers = set()
    for i in range(len(seq) - k + 1):
        parca = seq[i:i+k]
        kmers.add(parca)
    return kmers

def calc_distance(seq1, seq2):
    k1 = get_kmers(seq1)
    k2 = get_kmers(seq2)
    
    kesisim_sayisi = len(k1.intersection(k2))
    birlesim_sayisi = len(k1.union(k2))
    
    if birlesim_sayisi == 0:
        return 1
    else:
        uzaklik = 1 - (kesisim_sayisi / birlesim_sayisi)
        return uzaklik

def needleman_wunsch(seq1, seq2, gp=-1, m=1, mm=-1):
    n = len(seq1)
    l = len(seq2)
    
    sc = []
    for _ in range(n + 1):
        satir = [0] * (l + 1)
        sc.append(satir)
        
    for i in range(n + 1): 
        sc[i][0] = i * gp
    for j in range(l + 1): 
        sc[0][j] = j * gp

    for i in range(1, n + 1):
        for j in range(1, l + 1):
            if seq1[i-1] == seq2[j-1]:
                eslesme_puani = m
            else:
                eslesme_puani = mm
                
            match = sc[i-1][j-1] + eslesme_puani
            delete = sc[i-1][j] + gp
            insert = sc[i][j-1] + gp
            
            sc[i][j] = max(match, delete, insert)

    a1 = ""
    a2 = ""
    i = n
    j = l
    
    while i > 0 and j > 0:
        c = sc[i][j]
        
        if seq1[i-1] == seq2[j-1]:
            eslesme_puani = m
        else:
            eslesme_puani = mm
            
        if c == sc[i-1][j-1] + eslesme_puani:
            a1 = seq1[i-1] + a1
            a2 = seq2[j-1] + a2
            i -= 1
            j -= 1
        elif c == sc[i-1][j] + gp:
            a1 = seq1[i-1] + a1
            a2 = "-" + a2
            i -= 1
        else:
            a1 = "-" + a1
            a2 = seq2[j-1] + a2
            j -= 1
            
    while i > 0: 
        a1 = seq1[i-1] + a1
        a2 = "-" + a2
        i -= 1
    while j > 0: 
        a1 = "-" + a1
        a2 = seq2[j-1] + a2
        j -= 1
        
    return a1, a2

def get_consensus(a1, a2):
    consensus = ""
    for i in range(len(a1)):
        harf1 = a1[i]
        harf2 = a2[i]
        
        if harf1 == harf2:
            consensus += harf1
        elif harf1 == '-':
            consensus += harf2
        elif harf2 == '-':
            consensus += harf1
        else:
            consensus += harf1
            
    return consensus

def my_muscle_msa(seqs):
    if len(seqs) < 2: 
        return seqs
        
    s = seqs.copy()
    md = float('inf')
    pair = (0, 1)
    
    for i in range(len(s)):
        for j in range(i+1, len(s)):
            d = calc_distance(s[i], s[j])
            if d < md: 
                md = d
                pair = (i, j)

    a = s.pop(max(pair))
    b = s.pop(min(pair))
    
    hizali_a, hizali_b = needleman_wunsch(a, b)
    grp = [hizali_a, hizali_b]
    cons = get_consensus(grp[0], grp[1])

    while len(s) > 0:
        nxt = s.pop(0)
        ac, an = needleman_wunsch(cons, nxt)
        
        ng = []
        for _ in range(len(grp)):
            ng.append("")
            
        idx = 0
        for char in ac:
            if char == '-':
                for i in range(len(grp)): 
                    ng[i] += '-'
            else:
                for i in range(len(grp)): 
                    ng[i] += grp[i][idx]
                idx += 1
                
        ng.append(an)
        grp = ng
        cons = get_consensus(ac, an)
        
    return grp