//
// Created by zrobertson on 8/8/17.
//

#ifndef ECHELLESIMULATOR_TELESCOPE_H
#define ECHELLESIMULATOR_TELESCOPE_H


class Telescope {
public:
    Telescope();
    Telescope(double d, double f);
    ~Telescope();
    double get_area();
    double get_diameter();
    double get_focal_ratio();
private:
    double d;
    double f;

};


#endif //ECHELLESIMULATOR_TELESCOPE_H
