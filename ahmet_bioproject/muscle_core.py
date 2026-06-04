def get_kmers(seq, k=2):
    return set(seq[i:i+k] for i in range(len(seq)-k+1))

def calc_distance(seq1, seq2):
    k1, k2 = get_kmers(seq1), get_kmers(seq2)
    u = len(k1.union(k2))
    return 1 - (len(k1.intersection(k2)) / u) if u > 0 else 1

def needleman_wunsch(seq1, seq2, gp=-1, m=1, mm=-1):
    n, l = len(seq1), len(seq2)
    sc = [[0]*(l+1) for _ in range(n+1)]
    for i in range(n+1): sc[i][0] = i*gp
    for j in range(l+1): sc[0][j] = j*gp

    for i in range(1, n+1):
        for j in range(1, l+1):
            match = sc[i-1][j-1] + (m if seq1[i-1] == seq2[j-1] else mm)
            sc[i][j] = max(match, sc[i-1][j] + gp, sc[i][j-1] + gp)

    a1, a2, i, j = "", "", n, l
    while i > 0 and j > 0:
        c = sc[i][j]
        if c == sc[i-1][j-1] + (m if seq1[i-1] == seq2[j-1] else mm):
            a1 = seq1[i-1] + a1; a2 = seq2[j-1] + a2; i -= 1; j -= 1
        elif c == sc[i-1][j] + gp:
            a1 = seq1[i-1] + a1; a2 = "-" + a2; i -= 1
        else:
            a1 = "-" + a1; a2 = seq2[j-1] + a2; j -= 1
            
    while i > 0: a1 = seq1[i-1] + a1; a2 = "-" + a2; i -= 1
    while j > 0: a1 = "-" + a1; a2 = seq2[j-1] + a2; j -= 1
    return a1, a2

def get_consensus(a1, a2):
    return "".join(x if x == y else (y if x == '-' else (x if y == '-' else x)) for x, y in zip(a1, a2))

def my_muscle_msa(seqs):
    if len(seqs) < 2: return seqs
    s = seqs.copy()
    md, pair = float('inf'), (0, 1)
    for i in range(len(s)):
        for j in range(i+1, len(s)):
            d = calc_distance(s[i], s[j])
            if d < md: md, pair = d, (i, j)

    a, b = s.pop(max(pair)), s.pop(min(pair))
    grp = list(needleman_wunsch(a, b))
    cons = get_consensus(grp[0], grp[1])

    while s:
        nxt = s.pop(0)
        ac, an = needleman_wunsch(cons, nxt)
        ng, idx = [""] * len(grp), 0
        for char in ac:
            if char == '-':
                for i in range(len(grp)): ng[i] += '-'
            else:
                for i in range(len(grp)): ng[i] += grp[i][idx]
                idx += 1
        ng.append(an)
        grp, cons = ng, get_consensus(ac, an)
    return grp