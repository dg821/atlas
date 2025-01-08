#include "Dopri87.h"
#include <stdexcept>
#include <cmath>

Dopri87::Dopri87() : NumericalIntegrator(StepType::VariableStep) {
    initializeCoefficients();
}

void Dopri87::initializeCoefficients() {
    c << 0.0, 1.0/18.0, 1.0/12.0, 1.0/8.0, 5.0/16.0, 3.0/8.0, 59.0/400.0, 93.0/200.0,
         5490023248.0/9719169821.0, 13.0/20.0, 1201146811.0/1299019798.0, 1.0, 1.0;

    a.setZero();
    a(1, 0) = 1.0 / 18.0;
    a(2, 0) = 1.0 / 48.0;
    a(2, 1) = 1.0 / 16.0;
    a(3, 0) = 1.0 / 32.0;
    a(3, 2) = 3.0 / 32.0;
    a(4, 0) = 5.0 / 16.0;
    a(4, 2) = -75.0 / 64.0;
    a(4, 3) = 75.0 / 64.0;
    a(5, 0) = 3.0 / 80.0;
    a(5, 3) = 3.0 / 16.0;
    a(5, 4) = 3.0 / 20.0;
    a(6, 0) = 29443841.0 / 614563906.0;
    a(6, 3) = 77736538.0 / 692538347.0;
    a(6, 4) = -28693883.0 / 1125000000.0;
    a(6, 5) = 23124283.0 / 1800000000.0;
    a(7, 0) = 16016141.0 / 946692911.0;
    a(7, 3) = 61564180.0 / 158732637.0;
    a(7, 4) = 22789713.0 / 633445777.0;
    a(7, 5) = 545815736.0 / 2771057229.0;
    a(7, 6) = -180193667.0 / 1043307555.0;
    a(8, 0) = 39632708.0 / 573591083.0;
    a(8, 3) = -433636366.0 / 683701615.0;
    a(8, 4) = -421739975.0 / 2616292301.0;
    a(8, 5) = 100302831.0 / 723423059.0;
    a(8, 6) = 790204164.0 / 839813087.0;
    a(8, 7) = 800635310.0 / 3783071287.0;
    a(9, 0) = 246121993.0 / 1340847787.0;
    a(9, 3) = -37695042795.0 / 15268766246.0;
    a(9, 4) = -309121744.0 / 1061227803.0;
    a(9, 5) = -12992083.0 / 490766935.0;
    a(9, 6) = 6005943493.0 / 2108947869.0;
    a(9, 7) = 393006217.0 / 1396673457.0;
    a(9, 8) = 123872331.0 / 1001029789.0;
    a(10, 0) = -1028468189.0 / 846180014.0;
    a(10, 3) = 8478235783.0 / 508512852.0;
    a(10, 4) = 1311729495.0 / 1432422823.0;
    a(10, 5) = -10304129995.0 / 1701304382.0;
    a(10, 6) = -48777925059.0 / 3047939560.0;
    a(10, 7) = 15336726248.0 / 1032824649.0;
    a(10, 8) = -45442868181.0 / 3398467696.0;
    a(10, 9) = 3065993473.0 / 597172653.0;
    a(11, 0) = 185892177.0 / 718116043.0;
    a(11, 3) = -3185094517.0 / 667107341.0;
    a(11, 4) = -477755414.0 / 1098053517.0;
    a(11, 5) = -703635378.0 / 230739211.0;
    a(11, 6) = 5731566787.0 / 1027545527.0;
    a(11, 7) = 5232866602.0 / 850066563.0;
    a(11, 8) = -4093664535.0 / 808688257.0;
    a(11, 9) = 3962137247.0 / 1805957418.0;
    a(11, 10) = 65686358.0 / 487910083.0;
    a(12, 0) = 403863854.0 / 491063109.0;
    a(12, 3) = -5068492393.0 / 434740067.0;
    a(12, 4) = -411421997.0 / 543043805.0;
    a(12, 5) = 652783627.0 / 914296604.0;
    a(12, 6) = 11173962825.0 / 925320556.0;
    a(12, 7) = -13158990841.0 / 6184727034.0;
    a(12, 8) = 3936647629.0 / 1978049680.0;
    a(12, 9) = -160528059.0 / 685178525.0;
    a(12, 10) = 248638103.0 / 1413531060.0;

    b7 << 13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0, -808719846.0/976000145.0,
          1757004468.0/5645159321.0, 656045339.0/265891186.0, -3867574721.0/1518517206.0,
          465885868.0/322736535.0, 53011238.0/667516719.0, 2.0/45.0, 0.0;

    b8 << 14005451.0/335480064.0, 0.0, 0.0, 0.0, 0.0, -59238493.0/1068277825.0,
          181606767.0/758867731.0, 561292985.0/797845732.0, -1041891430.0/1371343529.0,
          760417239.0/1151165299.0, 118820643.0/751138087.0, -528747749.0/2220607170.0, 1.0/4.0;
}

Dopri87::StateVector Dopri87::step(const RightHandSide& f, double t, const StateVector& y, double h) const {
    auto [y_7, y_8] = computeStage(f, t, y, h);
    return y_8;  // Return the 8th order solution
}

std::vector<std::pair<double, Dopri87::StateVector>> Dopri87::integrate(
    const RightHandSide& f, double t_start, double t_end, const StateVector& y_start,
    double h, double h_min, double h_max) const {
    std::vector<std::pair<double, StateVector>> solution;
    solution.reserve(static_cast<size_t>((t_end - t_start) / h) + 1);
    solution.emplace_back(t_start, y_start);

    double t = t_start;
    StateVector y = y_start;

    while (t < t_end) {
        if (t + h > t_end) h = t_end - t;

        auto [y_new, h_new] = adaptiveStep(f, t, y, h, h_min, h_max);

        t += h;
        y = std::move(y_new);
        h = h_new;

        solution.emplace_back(t, y);
    }

    return solution;
}

std::pair<Dopri87::StateVector, Dopri87::StateVector> Dopri87::computeStage(
    const RightHandSide& f, double t, const StateVector& y, double h) const {
    std::array<Eigen::VectorXd, stages> k;
    k[0] = f(t, y);

    auto compute_k = [&](int i) {
        StateVector y_temp = y;
        for (int j = 0; j < i; ++j) {
            y_temp += h * a(i, j) * k[j];
        }
        return f(t + c(i) * h, y_temp);
    };

    for (int i = 1; i < stages; ++i) {
        k[i] = compute_k(i);
    }

    // Create a matrix from all k vectors
    Eigen::Matrix<double, Eigen::Dynamic, stages> K(y.size(), stages);
    for (int i = 0; i < stages; ++i) {
        K.col(i) = k[i];
    }

    StateVector y_7 = y + h * (K * b7);
    StateVector y_8 = y + h * (K * b8);

    return {y_7, y_8};
}

double Dopri87::estimateError(const StateVector& y_7, const StateVector& y_8) const {
    return ((y_8 - y_7).array() / (atol + rtol * y_8.array().abs().max(y_7.array().abs()))).matrix().norm() / std::sqrt(y_7.size());
}

std::pair<Dopri87::StateVector, double> Dopri87::adaptiveStep(
    const RightHandSide& f, double t, const StateVector& y, double h, double h_min, double h_max) const {
    const double safety = 0.9;
    const double minReduction = 0.1;
    const double maxGrowth = 5.0;

    while (true) {
        auto [y_7, y_8] = computeStage(f, t, y, h);
        double err = estimateError(y_7, y_8);

        double h_new = h * std::min(maxGrowth, std::max(minReduction, safety * std::pow(1.0 / err, 1.0 / 8.0)));
        h_new = std::min(h_max, std::max(h_min, h_new));

        if (err <= 1.0) {
            return {y_8, h_new};
        }

        h = h_new;
        if (h < h_min) {
            throw std::runtime_error("Step size became too small");
        }
    }
}

int Dopri87::getOrder() const {
    return order;
}

std::string Dopri87::getName() const {
    return "Dormand-Prince 8(7)";
}