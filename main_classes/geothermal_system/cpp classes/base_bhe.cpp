#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cmath>




class BaseBHE {
private:
    std::vector<ThermodynamicPoint> ideal_points;
    std::vector<ThermodynamicPoint> real_points;
    ThermodynamicPoint tmp_point;
    BHEGeometry geom;
    ReservoirProperties res_prop;
    std::string user_unit_system;
    // Placeholder for integrator, assuming its relevant methods and properties
    // IsobaricIntegral integrator;

public:
    // Constructor
    BaseBHE(const ThermodynamicPoint& input_point, const BHEGeometry& geometry = BHEGeometry(), const ReservoirProperties& reservoir_properties = ReservoirProperties())
    : geom(geometry), res_prop(reservoir_properties), user_unit_system(input_point.RPHandler.unit_system) {
        // Assuming input_point can be directly used
        ideal_points.push_back(input_point);
        for (int i = 0; i < 3; ++i) {
            ideal_points.push_back(input_point.duplicate());
        }
        real_points.push_back(input_point.duplicate());
        real_points.push_back(input_point.duplicate());

        // Assuming integrator setup as needed
        // integrator = IsobaricIntegral(ideal_points);
    }

    // Example property methods
    double ideal_exergy_efficiency() const {
        double dh = ideal_points.back().get_variable("H") - ideal_points.front().get_variable("H");
        double ds = ideal_points.back().get_variable("S") - ideal_points.front().get_variable("S");
        double t_0 = ideal_points.front().get_variable("T");
        double t_rocks = t_0 + res_prop.grad * geom.depth;
        double dex_rocks = dh * (1 - t_0 / t_rocks);
        return (dh - t_0 * ds) / dex_rocks;
    }

    ThermodynamicPoint get_input_point() const {
        return ideal_points.front().get_alternative_unit_system(user_unit_system);
    }

    void set_input_point(const ThermodynamicPoint& input_point) {
        ideal_points.front() = input_point.get_alternative_unit_system("MASS BASE SI");
    }

    std::vector<ThermodynamicPoint> get_ideal_points() const {
        std::vector<ThermodynamicPoint> result;
        std::transform(ideal_points.begin(), ideal_points.end(), std::back_inserter(result),
                       [this](const ThermodynamicPoint& point) { return point.get_alternative_unit_system(user_unit_system); });
        return result;
    }

    std::vector<ThermodynamicPoint> get_real_points() const {
        std::vector<ThermodynamicPoint> result;
        result.push_back(ideal_points.front().get_alternative_unit_system(user_unit_system));
        result.push_back(ideal_points.at(1).get_alternative_unit_system(user_unit_system));

        for (const auto& point : real_points) {
            result.push_back(point.get_alternative_unit_system(user_unit_system));
        }
        return result;
    }

    // Additional methods would need to be similarly translated
};
