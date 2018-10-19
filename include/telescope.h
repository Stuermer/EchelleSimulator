#ifndef ECHELLESIMULATOR_TELESCOPE_H
#define ECHELLESIMULATOR_TELESCOPE_H

/*!
 * \class Telescope
 * \brief Telescope class
 */
class Telescope {
public:
    /**
     * Default constructor. Creates telescope with a 1 m^2 mirror
     */
    Telescope();

    /**
     * Constructor
     * @param d primary mirror diameter [m]
     * @param f total focal length [m]
     */
    Telescope(double d, double f);

    /**
     * Deconstructor
     */
    ~Telescope() {};

    ///@return effective light collecting area [m^2]
    double get_area();

    /// @return mirror diameter
    double get_diameter();

    /// @return focal ratio
    double get_focal_ratio();

private:
    /// mirror diameter [m]
    double d;

    /// total focal length [m]
    double f;

};


#endif //ECHELLESIMULATOR_TELESCOPE_H
