Point(1336) = {-0.0004018925337573828, 0.05177306826503763, 0, 5e-05};
// Import a line only that case when the start and the end point is still xists
Point(1) = {0, 0, 0, lc}; // Line(2) = {2,3}; is also an invalid comment
Point(2) = {.1, 0,  0, lc};
Line(1) = {1, 2}; // this should be in the results
// invalid line object, commented out
// Line(2) = {2,3};
//Circle(60) = {33, 291, 57};
