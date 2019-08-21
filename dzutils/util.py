def read_file_lines(filename):
    with open(filename, "r") as myfile:
        lines = myfile.readlines()
    return lines


# TODO: DOES not yet handle comments after the options!
def read_flag_file(filename):
    """Reads the flag file, ignoring comments"""
    lines = read_file_lines(filename)
    # filter the lines
    lines = [l for l in lines if l.startswith("-")]
    return " ".join(lines)


def pythonify(json_data):

    correctedDict = {}

    for key, value in json_data.items():
        if isinstance(value, list):
            value = [
                pythonify(item) if isinstance(item, dict) else item
                for item in value
            ]
        elif isinstance(value, dict):
            value = pythonify(value)
        try:
            key = int(key)
        except Exception as ex:
            pass
        correctedDict[key] = value

    return correctedDict
