Point(1336) = {-0.0004018925337573828, 0.05177306826503763, 0, 5e-05};
// Import a line only that case when the start and the end point is still xists
Point(1) = {0, 0, 0, lc};
Point(2) = {.1, 0,  0, lc}; Point(3) = {.1, 0.1,  0, lc};
Line(1) = {1, 2};
Circle(60) = {1, 2, 1336};
Point(101) = {2.5, 36, 0.0, lc};
Point(102) = {1.3, 35.6, 0.0, lc};
Point(103) = {1.3, 34.4, 0.0, lc};
Circle(160) = {101, 102, 103};
