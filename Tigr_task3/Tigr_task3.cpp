#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

struct InputData {
    double p;
    int n;
    double rho;
    int T;
    double delta;
    double eta;
    vector<double> eta_vec;
    double epsilon;
    double alpha;
    double beta;
    double gamma;
    vector<vector<double>> c0;
    vector<vector<double>> pi;
    vector<vector<int>> A;
};

InputData readInput() {
    InputData D;
    if (!(cin >> D.p)) { cerr << "Input error\n"; exit(1); }
    cin >> D.n >> D.rho >> D.T >> D.delta >> D.eta;

    D.eta_vec.assign(D.n, 0.0);
    for (int i = 0; i < D.n; ++i) cin >> D.eta_vec[i];

    cin >> D.epsilon >> D.alpha >> D.beta >> D.gamma;

    D.c0.assign(D.n, vector<double>(D.T + 1, 0.0));
    for (int i = 0; i < D.n; ++i) cin >> D.c0[i][0];

    D.pi.assign(D.n, vector<double>(D.n, 0.0));
    for (int i = 0; i < D.n; ++i)
        for (int j = 0; j < D.n; ++j)
            cin >> D.pi[i][j];

    D.A.assign(D.n, vector<int>(D.n, 0));
    return D;
}

void recalcA_dynamic(InputData& D, const vector<double>& Y_current, const vector<double>& Phi_T, int T_val) {
    int n = D.n;
    vector<vector<int>> newA(n, vector<int>(n, 0));

    double rho_power = pow(D.rho, T_val - 1);
    if (rho_power < 1e-15) rho_power = 1e-15;
    double inv_rho_power = -1.0 / rho_power;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double f = inv_rho_power * (D.beta - D.gamma) * Y_current[j] * Phi_T[i];
            if (D.pi[i][j] < f) {
                newA[i][j] = 1;
            }
        }
    }
    D.A = move(newA);
}

int main() {
    setlocale(LC_ALL, "Russian");
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cerr << "Введите данные:\n";
    InputData D = readInput();
    int n = D.n;
    int T = D.T;

    vector<vector<double>> C = D.c0;
    vector<vector<double>> U(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> Phi(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> Y(n, vector<double>(T + 1, 0.0));
    vector<double> total_profit(n, 0.0);

    // Терминальное условие Phi[i][T]
    vector<double> Phi_T(n);
    for (int i = 0; i < n; ++i) {
        Phi_T[i] = -pow(D.rho, T) * D.eta;
        Phi[i][T] = Phi_T[i];
    }

    // Основной цикл
    for (int t = 0; t < T; ++t) {
        // 1. U
        for (int i = 0; i < n; ++i) {
            double sum_other = 0.0;
            for (int j = 0; j < n; ++j) if (j != i) sum_other += C[j][t];
            U[i][t] = (D.p - (n + 1) * C[i][t] + sum_other) / (n + 1);
            if (U[i][t] < 0.0) U[i][t] = 0.0;
        }

        // 2. Phi[t]
        for (int i = 0; i < n; ++i) {
            double rho_t = pow(D.rho, t);
            Phi[i][t] = -rho_t * U[i][t] + D.delta * Phi[i][t + 1];
        }

        // 3. Y[t]
        for (int i = 0; i < n; ++i) {
            double rho_t = pow(D.rho, t);
            if (rho_t * D.epsilon > 1e-12) {
                Y[i][t] = -D.alpha * Phi[i][t] / (rho_t * D.epsilon);
            }
            else {
                Y[i][t] = 0.0;
            }
            if (Y[i][t] < 0.0) Y[i][t] = 0.0;
        }

        // 4. Обновляем A на шаге t, используя:
        //    - Y_current = Y[t]
        //    - Phi_T = Phi[i][T] (фиксировано!)
        vector<double> Y_current(n);
        for (int i = 0; i < n; ++i) Y_current[i] = Y[i][t];
        recalcA_dynamic(D, Y_current, Phi_T, T);

        // Вывод A
        cout << "Матрица A на t=" << t << ":\n";
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j > 0) cout << "\t";
                cout << D.A[i][j];
            }
            cout << "\n";
        }
        cout << "---------------------------\n";

        // 5. Обновление C
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                int g_ij = D.A[i][j];
                int g_ji = D.A[j][i];
                double interaction = (D.beta * g_ij * g_ji + D.gamma * (1 - g_ij * g_ji));
                sum += interaction * Y[j][t];
            }
            double next_c = D.delta * C[i][t] - D.alpha * Y[i][t] - sum;
            if (next_c < 0.0) next_c = 0.0;
            C[i][t + 1] = next_c;
        }
    }

    for (int i = 0; i < n; ++i) Y[i][T] = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum_other = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) sum_other += C[j][T];
        U[i][T] = (D.p - (n + 1) * C[i][T] + sum_other) / (n + 1);
        if (U[i][T] < 0.0) U[i][T] = 0.0;
    }

    for (int t = 0; t <= T; ++t) {
        double Q = 0.0;
        for (int j = 0; j < n; ++j) Q += U[j][t];
        double P_t = D.p - Q;
        for (int i = 0; i < n; ++i) {
            double profit_t = P_t * U[i][t] - 0.5 * D.epsilon * U[i][t] * U[i][t];
            total_profit[i] += profit_t;
        }
    }

    cout.setf(std::ios::fixed);
    cout << setprecision(6);

    for (int i = 0; i < n; ++i) {
        cout << "===========================\n";
        cout << "Фирма " << (i + 1) << "\n";
        cout << "t:\t";
        for (int t = 0; t <= T; ++t) { if (t > 0) cout << "\t"; cout << t; }
        cout << "\n";

        cout << "u_" << (i + 1) << "(t)\t";
        for (int t = 0; t <= T; ++t) { if (t > 0) cout << "\t"; cout << U[i][t]; }
        cout << "\n";

        cout << "y_" << (i + 1) << "(t)\t";
        for (int t = 0; t <= T; ++t) { if (t > 0) cout << "\t"; cout << Y[i][t]; }
        cout << "\n";

        cout << "прибыль в игре:\t" << total_profit[i] << "\n";
        cout << "===========================\n\n";
    }

    cout << "Свод по периодам (t, P_t, Q_t):\n";
    for (int t = 0; t <= T; ++t) {
        double Q = 0.0;
        for (int j = 0; j < n; ++j) Q += U[j][t];
        double P_t = D.p - Q;
        cout << "t=" << t << "  P=" << P_t << "  Q=" << Q << "\n";
    }

    return 0;
}