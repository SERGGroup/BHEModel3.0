#ifndef RESERVOIRPROPERTIES_H
#define RESERVOIRPROPERTIES_H

class ReservoirProperties {
public:
    double grad;
    double evaluate_rel_resistance(double time, double d) const; // Assume this is a method
    // Assume more properties and methods as required
};

#endif