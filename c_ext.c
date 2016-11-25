int str_compare(int len1, char* seq1, int* qual1,
                int len2, char* seq2, int* qual2,
                int qual) {
    int i;

    if (len1 != len2)
        return 0;

    if (*seq1 != *seq2) // rule of the origin code; I dont know why
        return 0;

    for (i = 1; i < len1; i++) {
        if (seq1[i] == seq2[i] || qual1[i] < qual || qual2[i] < qual) {
            continue;
        }
        return 0;
    }

    return 1;
}
