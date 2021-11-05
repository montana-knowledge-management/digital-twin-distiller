from math import log, pow


def inductancefromimages(xi, yi, xw1, yw1, xw2, yw2):
    constant = 0.2

    d1 = pow(pow(xi - xw1, 2) + pow(yi - yw1, 2), 0.5)
    d2 = pow(pow(xi - xw2, 2) + pow(yi - yw2, 2), 0.5)

    return constant * log(d1 / d2)


def simpleCalc(x_1, y_1, x_2, y_2, r, layernum):

    w = 0.1
    h = 0.3
    result = 0

    for i in range(-layernum, layernum):
        for j in range(-layernum, layernum):
            result += inductancefromimages(x_2 + 2 * i * w, y_2 + 2 * j * h, 0, y_1, x_1 - r, y_1)
            result -= inductancefromimages(-x_2 + 2 * i * w, y_2 + 2 * j * h, 0, y_1, x_1 - r, y_1)
            result -= inductancefromimages(x_2 + 2 * i * w, -y_2 + 2 * j * h, 0, y_1, x_1 - r, y_1)
            result += inductancefromimages(-x_2 + 2 * i * w, -y_2 + 2 * j * h, 0, y_1, x_1 - r, y_1)

    return result


def main():
    r = 0.35 * 1e-3
    print("r", r)
    res = 0
    for i in range(27):
        for j in range(27):
            if i == j:
                res += (
                    simpleCalc(
                        1.0 * 1e-3 + r,
                        26.5 * 1e-3 - i * 1e-3,
                        1.0e-3 + r,
                        (26.5 * 1e-3) - j * 1e-3,
                        r,
                        100,
                    )
                    / 2.0
                )
                print(
                    i,
                    simpleCalc(
                        1.0 * 1e-3,
                        26.5 * 1e-3 - i * 1e-3,
                        1.0 * 1e-3 + r,
                        (26.5 * 1e-3) - j * 1e-3,
                        r,
                        100,
                    ),
                )

    print("Total:", res)


if __name__ == "__main__":
    main()
