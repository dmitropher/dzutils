def string_to_fasta(s, name=">FASTA Sequence", width=80):
    """
    Returns a FASTA formated string from the given input string (Must be AA or NA, no comments)

    optional:
    name - can specify a header for the fasta file, must begin with > and be shorter than the width specified
    width - character width to assign to the FASTA, 0 just writes the string to one line
    """
    if width == 0:
        return name + "\n" + s
    fasta = (
        name
        + "\n"
        + "\n".join([s[i : i + width] for i in range(0, len(s), width)])
    )
    return fasta
