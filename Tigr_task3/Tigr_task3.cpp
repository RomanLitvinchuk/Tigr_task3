#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

// Параметры сходимости и итераций
const double EPSILON_CONVERGENCE = 1e-6;
const int MAX_ITERATIONS = 20000;

struct Parameters {
    double p;             // Максимальная стоимость единицы товара
    int n;                // Количество фирм
    double rho;           // Коэф. дисконтирования [0; 1]
    int T;                // Терминальный момент времени
    double delta;         // Коэф. сохранения издержек
    double eta; // Ликвидационная стоимость (max рыночная стоимость)
    vector<double> eta_vec;
    double epsilon;       // Параметр (цена усилий)
    double alpha;         // Эффективность инвестиций
    double beta;          // Эффективность сотрудничества
    double gamma;         // Эффективность конкуренции (без сотрудничества)
    vector<double> c0;    // Начальные удельные издержки
    vector<vector<double>> pi; // Матрица издержек сотрудничества
};

struct GameState {
    // Вектора состояний размером [n][T+1]
    vector<vector<double>> c;   // Удельные издержки
    vector<vector<double>> u;   // Объем товара
    vector<vector<double>> y;   // Инвестиционные усилия
    vector<vector<double>> phi; // Сопряженная переменная
    vector<vector<double>> profit; // Прибыль
    vector<double> total_discounted_profit; // Суммарная дисконтированная прибыль

    vector<vector<vector<int>>> g;

    void resize(int n, int T) {
        c.assign(n, vector<double>(T + 1, 0.0));
        u.assign(n, vector<double>(T + 1, 0.0));
        y.assign(n, vector<double>(T + 1, 0.0));
        phi.assign(n, vector<double>(T + 1, 0.0));
        profit.assign(n, vector<double>(T + 1, 0.0));
        total_discounted_profit.assign(n, 0.0); 
        g.assign(T + 1, vector<vector<int>>(n, vector<int>(n, 0)));
    }
};

void solveGame(const Parameters& params, GameState& state) {
    int n = params.n;
    int T = params.T;

    state.resize(n, T);
    for (int i = 0; i < n; ++i) {
        state.c[i][0] = params.c0[i];
    }

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        double max_diff = 0.0;

        // 1. Расчет производства u(t) (t < T)
        for (int t = 0; t < T; ++t) {
            double sum_c = 0.0;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    sum_c += state.c[j][t];
                    //double val = (params.p - (n + 1) * state.c[i][t] + sum_c) / (n + 1);
                    //state.u[i][t] = max(0.0, val);
                }
                state.u[i][t] = (params.p - (n + 1) * state.c[i][t] + sum_c) / (n + 1);
                state.u[i][t] = max(0.0, state.u[i][t]);
            }
        }
        for (int i = 0; i < n; ++i) state.u[i][T] = 0.0;

        // 2. Расчет сопряженной переменной phi(t) 
        for (int i = 0; i < n; ++i) {
            state.phi[i][T] = -pow(params.rho, T) * params.eta;

            for (int t = T - 1; t >= 0; --t) {
                double next_phi = -pow(params.rho, t) * state.u[i][t] + params.delta * state.phi[i][t + 1];
                max_diff = max(max_diff, abs(next_phi - state.phi[i][t]));
                state.phi[i][t] = next_phi;
            }
        }

        // Расчет инвестиций y(t) и матрицы смежности g(t) (t < T)
        for (int t = 0; t < T; ++t) {
            for (int i = 0; i < n; ++i) {
                double rho_t = pow(params.rho, t);
                double denom = rho_t * params.epsilon;
                double val = 0.0;

                if (abs(denom) > 1e-9) {
                    val = (-params.alpha / denom) * state.phi[i][t + 1];
                }
                state.y[i][t] = max(0.0, val);
            }

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    double benefit = -state.phi[i][t + 1] * (params.beta - params.gamma) * state.y[j][t];

                    if (benefit > params.pi[i][j]) {
                        state.g[t][i][j] = 1;
                    }
                    else {
                        state.g[t][i][j] = 0;
                    }
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            state.y[i][T] = 0.0;
            for (int j = 0; j < n; ++j) state.g[T][i][j] = 0;
        }

        // 4. Расчет динамики издержек c(t+1)
        for (int t = 0; t < T; ++t) {
            for (int i = 0; i < n; ++i) {
                double sum_interaction = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    int cooperation = state.g[t][i][j] * state.g[t][j][i];
                    double term = (params.beta * cooperation + params.gamma * (1 - cooperation)) * state.y[j][t];
                    sum_interaction += term;
                }

                double next_c = params.delta * state.c[i][t] - params.alpha * state.y[i][t] - sum_interaction;
                next_c = max(0.0, next_c);

                max_diff = max(max_diff, abs(next_c - state.c[i][t + 1]));
                state.c[i][t + 1] = next_c;
            }
        }

        if (max_diff < EPSILON_CONVERGENCE) {
            break;
        }
    }

    // Расчет мгновенной и суммарной дисконтированной прибыли
    for (int i = 0; i < n; ++i) {
        state.total_discounted_profit[i] = 0.0; 
    }
    for (int t = 0; t < T; ++t) {
        double coop_wastes = 0.0;
        double total_uj = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) 
            {
                if (i != j) { coop_wastes += params.pi[i][j] * state.g[t][i][j] * state.g[t][j][i]; }
                total_uj += state.u[j][t];
            }
            state.total_discounted_profit[i] += pow(params.rho, t) * ((params.p - state.c[i][t] - total_uj) * state.u[i][t] - params.epsilon / 2 * pow(state.y[i][t], 2) - coop_wastes) + pow(params.rho, T) * (params.eta_vec[i] - params.eta * state.c[i][T]);
        }

    }
}

void printResults(const Parameters& params, const GameState& state) {
    cout << fixed << setprecision(4);

    for (int i = 0; i < params.n; ++i) {
        cout << "\n==========================================" << endl;
        cout << "Table for Firm " << (i + 1) << ":" << endl;
        cout << "------------------------------------------" << endl;

        cout << left << setw(15) << "t";
        for (int t = 0; t <= params.T; ++t) cout << setw(10) << t;
        cout << endl << "------------------------------------------" << endl;

        cout << left << setw(15) << ("u_" + to_string(i + 1) + "(t)");
        for (int t = 0; t <= params.T; ++t) cout << setw(10) << state.u[i][t];
        cout << endl;

        cout << left << setw(15) << ("y_" + to_string(i + 1) + "(t)");
        for (int t = 0; t <= params.T; ++t) cout << setw(10) << state.y[i][t];
        cout << endl << "------------------------------------------" << endl;

        cout << left << setw(15) << "Profit (flow)";
        for (int t = 0; t <= params.T; ++t) cout << setw(10) << state.profit[i][t];
        cout << endl;

        cout << "\n>>> Total Discounted Profit: " << state.total_discounted_profit[i] << endl;
    }

    cout << "\n\n=== Adjacency Matrices (Is cooperation profitable?) ===" << endl;
    for (int t = 0; t <= params.T; ++t) {
        cout << "Time t = " << t << (t == params.T ? " (Terminal state - all zeros)" : "") << ":" << endl;
        cout << "   ";
        for (int k = 0; k < params.n; ++k) cout << "F" << k + 1 << " ";
        cout << endl;

        for (int i = 0; i < params.n; ++i) {
            cout << "F" << i + 1 << " ";
            for (int j = 0; j < params.n; ++j) {
                cout << setw(2) << state.g[t][i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

int main() 
{
    Parameters params;

    cout << "Enter p (max price of product): ";
    cin >> params.p;

    cout << "Enter n (number of firms): ";
    cin >> params.n;

    cout << "Enter rho (discount factor [0; 1]): ";
    cin >> params.rho;

    cout << "Enter T (terminal time, integer): ";
    cin >> params.T;

    cout << "Enter delta (cost retention coeff [0; 1]): ";
    cin >> params.delta;

    cout << "Enter eta (liquidity coeff) for " << params.n << " firms:" << endl;
    cin >> params.eta;

    cout << "Enter eta_i for each firm:" << endl;
    params.eta_vec.resize(params.n);
    for (int i = 0; i < params.n; ++i)
    {
        cout << "  eta[" << i + 1 << "]:" << endl;
        cin >> params.eta_vec[i];
    }
    cout << "Enter epsilon (cost param), alpha (inv eff), beta (coop eff), gamma (non-coop eff): " << endl;
    cout << "Format: eps alpha beta gamma" << endl;
    cin >> params.epsilon >> params.alpha >> params.beta >> params.gamma;

    cout << "Enter initial costs c(0) for " << params.n << " firms:" << endl;

    params.c0.resize(params.n);
    for (int i = 0; i < params.n; ++i) {
        cout << "  c_" << i + 1 << "(0): ";
        cin >> params.c0[i];
    }

    cout << "Enter cooperation costs matrix pi (" << params.n << "x" << params.n << "):" << endl;
    cout << "(Row i Column j = cost for firm i to cooperate with j)" << endl;
    params.pi.assign(params.n, vector<double>(params.n));
    for (int i = 0; i < params.n; ++i) {
        for (int j = 0; j < params.n; ++j) {
            cout << "  pi[" << i + 1 << "][" << j + 1 << "]: ";
            cin >> params.pi[i][j];
        }
    }

    GameState state;
    solveGame(params, state);
    printResults(params, state);

    return 0;
}