for line in open("earth_spec.ini"):
    if line[0] == "#": continue
    elements = line.split()
    if elements[0] == "continue":
        print("continue {} {}".format(elements[1], float(elements[2]) - 4))
    elif elements[0] == "intensity":
        print("intensity {}".format(float(elements[1])))
    else:
        print(line)
