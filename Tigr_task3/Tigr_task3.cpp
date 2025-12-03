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

    return D;
}

// ---------------------- Генерация матрицы A по минимизации издержек ----------------------
void generateOptimalA(InputData& D) {
    int n = D.n;
    D.A.assign(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        int best_j = -1;
        double min_cost = 1e18;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double coopCost = D.pi[i][j]; // простая минимизация издержек
            double totalCost = D.c0[i][0] + coopCost;
            if (totalCost < min_cost) {
                min_cost = totalCost;
                best_j = j;
            }
        }
        if (best_j != -1) D.A[i][best_j] = 1;
    }

    cout << "Матрица минимизации издержек A:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j > 0) cout << "\t";
            cout << D.A[i][j];
        }
        cout << "\n";
    }
    cout << "-------------------------------------\n";
}

// ---------------------- Основная программа ----------------------
int main() {
    setlocale(LC_ALL, "Russian");
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cerr << "Введите данные (без матрицы A, она будет сформирована автоматически):\n";
    InputData D = readInput();
    int n = D.n;
    int T = D.T;

    generateOptimalA(D);

    vector<vector<double>> C = D.c0;
    vector<vector<double>> U(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> Phi(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> Y(n, vector<double>(T + 1, 0.0));
    vector<double> total_profit(n, 0.0);

    // ---------------------- Расчет u_i(t) ----------------------
    for (int t = 0; t <= T; ++t) {
        for (int i = 0; i < n; ++i) {
            double sum_other = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) sum_other += C[j][t];
            }
            U[i][t] = (D.p - (n + 1) * C[i][t] + sum_other) / (n + 1);
            if (U[i][t] < 0.0) U[i][t] = 0.0; // ставим 0 если отрицательное
        }
    }

    // ---------------------- Расчет phi_i_i(t) ----------------------
    for (int i = 0; i < n; ++i) Phi[i][T] = -pow(D.rho, T) * D.eta_vec[i];
    for (int t = T - 1; t >= 0; --t) {
        double rho_t = pow(D.rho, t);
        for (int i = 0; i < n; ++i) {
            Phi[i][t] = -rho_t * U[i][t] + D.delta * Phi[i][t + 1];
        }
    }

    // ---------------------- Расчет y_i(t) ----------------------
    for (int t = 0; t < T; ++t) {
        double rho_t = pow(D.rho, t);
        for (int i = 0; i < n; ++i) {
            Y[i][t] = D.alpha * Phi[i][t] / (rho_t * D.epsilon);
            double nextc = C[i][t] - D.alpha * U[i][t];
            if (nextc < 0.0) nextc = 0.0;
            C[i][t + 1] = nextc;
        }
    }
    for (int i = 0; i < n; ++i) Y[i][T] = 0.0;

    // ---------------------- Пересчет u_i(t) после обновления C ----------------------
    for (int t = 0; t <= T; ++t) {
        for (int i = 0; i < n; ++i) {
            double sum_other = 0.0;
            for (int j = 0; j < n; ++j) if (j != i) sum_other += C[j][t];
            U[i][t] = (D.p - (n + 1) * C[i][t] + sum_other) / (n + 1);
            if (U[i][t] < 0.0) U[i][t] = 0.0;
        }
    }

    // ---------------------- Расчет прибыли ----------------------
    for (int t = 0; t <= T; ++t) {
        double Q = 0.0;
        for (int j = 0; j < n; ++j) Q += U[j][t];
        double P_t = D.p - Q;
        for (int i = 0; i < n; ++i) {
            double profit_t = P_t * U[i][t] - 0.5 * D.epsilon * U[i][t] * U[i][t];
            total_profit[i] += profit_t;
        }
    }

    // ---------------------- Вывод ----------------------
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