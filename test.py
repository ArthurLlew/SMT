def is_point_in_D(x,y):
    if (y >= 0 and ((4/3)*x + 4 >= y) and ((-4/3)*x + 4 >= y)):
        return True
    return False


print(is_point_in_D(-2.7,0.2))


def get_segment_crossing(x1, y1, x2, y2, x3, y3, x4, y4):
    d = (x4 - x3)*(y2 - y1) - (x2 - x1)*(y4 - y3)

    s0 = ((x4 - x3)*(y3 - y1) - (x3 - x1)*(y4 - y3))/d
    t0 = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))/d

    if (0 < s0) and (0 < t0) and (s0 < 1) and (t0 < 1):
        x, y = x1 + s0*(x2 - x1), y1 + s0*(y2 - y1)
        return x, y
    else:
        return False

print(get_segment_crossing(-3, 0, 0, 4, -2.7, 0.6, -2.7, 0.2))


def get_segment_crossing_with_D(x1, y1, x2, y2):
    crossing = get_segment_crossing(-3, 0, 0, 4, x1, y1, x2, y2)
    if (crossing != False):
        return crossing
    crossing = get_segment_crossing(0, 4, 3, 0, x1, y1, x2, y2)
    if (crossing != False):
        return crossing
    if ((y1 == 0) and (y2 == 0)):
        if (x1 < -3) or (x2 < -3):
            return -3, 0
        if (x1 > 3) or (x2 > 3):
            return 3, 0
    crossing = get_segment_crossing(-3, 0, 3, 0, x1, y1, x2, y2)
    if (crossing != False):
        return crossing
    return False

print(get_segment_crossing_with_D(-2.7, 0.6, -2.7, 0.2))

def get_polygon_area(points: list):
    points.append(points[0])

    a = 0
    for i in range(len(points) - 1):
        a += points[i][1] * points[i + 1][0]

    b = 0
    for i in range(len(points) - 1):
        b += points[i][0] * points[i + 1][1]

    return abs(a - b)/2

points = []
points.append((-2.625000, 0.500000))
points.append((-2.550000, 0.500000))
points.append((-2.550000, 0.300000))
points.append((-2.775000, 0.300000))
print(get_polygon_area(points))