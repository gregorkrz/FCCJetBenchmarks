# Use these functions to format numbers in a root df column-friendly name


def point_format(number):
    return str(number).replace(".", "p")


def neg_format(number):
    # put n5 for -5
    # put n5p5 for -5.5
    if number < 0:
        return point_format("n{}".format(abs(number)))
    else:
        return point_format(number)
