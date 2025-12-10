#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>


using namespace std;

// Портирование C# Solver -> C++
// Структуры данных, имена и логика максимально сохранены.

struct InputData {
    int num;                // number of firms
    int time;               // time horizon (T)
    double discount;        // discount factor (rho)
    double liquidity;       // liquidity (eta in previous code)
    vector<vector<double>> unitCosts; // unitCosts[num][time+1] (we'll use [:][0] as initial)
    double personalInvestments; // personalInvestments (alpha in previous notation)
    double conversion;      // conversion (epsilon in previous notation)
    double maxProductCost;  // maxProductCost (p)
    vector<vector<double>> cooperationCosts; // cooperationCosts[num][num] (pi)
    double friendsInvestments; // friendsInvestments (beta)
    double notFriendsInvestments; // notFriendsInvestments (gamma)
    vector<double> maxCompanyValues; // maxCompanyValues[num] (eta_vec)
    double aging;           // aging (delta)
};

struct OutputData {
    vector<vector<double>> y; // [num][time+1]
    vector<vector<double>> u; // [num][time+1]
    int time;
    vector<double> profits; // [num]
    vector<vector<vector<int>>> g; // [num][num][time+1]
    int num;
};

OutputData Solve(const InputData& data) {
    int n = data.num;
    int T = data.time;

    // Allocate arrays
    vector<vector<double>> c(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> phi(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> y(n, vector<double>(T + 1, 0.0));
    vector<vector<double>> u(n, vector<double>(T + 1, 0.0));
    vector<vector<vector<int>>> g(n, vector<vector<int>>(n, vector<int>(T + 1, 0)));

    // 1) initial c: c[i,t] = c[i,0] for all t (как в C#)
    for (int i = 0; i < n; ++i) {
        for (int t = 0; t <= T; ++t) {
            c[i][t] = data.unitCosts[i][0];
        }
    }

    // 2) phi[T] = - discount^T * liquidity
    for (int i = 0; i < n; ++i) {
        phi[i][T] = -pow(data.discount, T) * data.liquidity;
    }

    // 3) y at time T-1 initial (as in C#)
    if (T - 1 >= 0) {
        for (int i = 0; i < n; ++i) {
            double denom = pow(data.discount, T - 1) * data.conversion;
            if (fabs(denom) < 1e-300) y[i][T - 1] = 0.0;
            else y[i][T - 1] = -data.personalInvestments / denom * phi[i][T];
        }
    }

    // 4) initial u: u[i,0] based on unitCosts[:,0], and set u[:,t]=u[:,0] for t>0 (as in C# init)
    for (int i = 0; i < n; ++i) {
        double sumOtherCosts = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) sumOtherCosts += data.unitCosts[j][0];
        double val = (data.maxProductCost - (n + 1) * data.unitCosts[i][0] + sumOtherCosts) / (n + 1);
        u[i][0] = val;
        for (int t = 1; t <= T; ++t) u[i][t] = u[i][0];
    }

    // 5) initial g at time T-1 (as in C#)
    if (T - 1 >= 0) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double rhs = 0.0;
                double denom = pow(data.discount, T - 1);
                if (fabs(denom) > 1e-300) {
                    rhs = -1.0 / denom * (data.friendsInvestments - data.notFriendsInvestments) * y[j][T - 1] * phi[i][T];
                }
                g[i][j][T - 1] = (data.cooperationCosts[i][j] < rhs) ? 1 : 0;
            }
        }
    }

    // Iterative process parameters
    const double eps = 1e-100; // as in C# code
    const int maxIterations = 1000;

    // Containers to store old values for convergence check
    vector<vector<double>> old_c(n, vector<double>(T + 1));
    vector<vector<double>> old_phi(n, vector<double>(T + 1));
    vector<vector<double>> old_y(n, vector<double>(T + 1));
    vector<vector<double>> old_u(n, vector<double>(T + 1));

    for (int iter = 0; iter < maxIterations; ++iter) {
        // copy current to old
        old_c = c;
        old_phi = phi;
        old_y = y;
        old_u = u;

        // --------- Backward sweep ----------
        for (int t = T - 1; t >= 0; --t) {
            for (int i = 0; i < n; ++i) {
                double denom = pow(data.discount, t) * data.conversion;
                if (fabs(denom) < 1e-300) y[i][t] = 0.0;
                else y[i][t] = -data.personalInvestments / denom * phi[i][t + 1];
            }

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (i == j) {
                        g[i][j][t] = 0;
                        continue;
                    }
                    double rhs = 0.0;
                    double denom = pow(data.discount, t);
                    if (fabs(denom) > 1e-300) {
                        rhs = -1.0 / denom * (data.friendsInvestments - data.notFriendsInvestments) * y[j][t] * phi[i][t + 1];
                    }
                    g[i][j][t] = (data.cooperationCosts[i][j] < rhs) ? 1 : 0;
                }
            }

            for (int i = 0; i < n; ++i) {
                phi[i][t] = -pow(data.discount, t) * u[i][t] + data.aging * phi[i][t + 1];
            }
        }

        // --------- Forward sweep ----------
        for (int t = 0; t <= T; ++t) {
            for (int i = 0; i < n; ++i) {
                // sumG as in C#
                double sumG = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    int gij = g[i][j][t];
                    int gji = g[j][i][t];
                    int coop_both = gij * gji;
                    double term = coop_both * y[j][t] * data.friendsInvestments
                        + (1 - coop_both) * y[j][t] * data.notFriendsInvestments;
                    sumG += term;
                }

                if (t < T) {
                    c[i][t + 1] = data.aging * c[i][t] - data.personalInvestments * y[i][t] - sumG;
                }
            }

            // recompute u[i,t] using c[:,t]
            double sumC = 0.0;
            for (int j = 0; j < n; ++j) sumC += c[j][t];
            for (int i = 0; i < n; ++i) {
                u[i][t] = (data.maxProductCost - (n + 1) * c[i][t] + sumC) / (n + 1);
            }
        }

        // Convergence check
        double maxChange = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int t = 0; t <= T; ++t) {
                maxChange = max(maxChange, fabs(c[i][t] - old_c[i][t]));
                maxChange = max(maxChange, fabs(phi[i][t] - old_phi[i][t]));
                maxChange = max(maxChange, fabs(y[i][t] - old_y[i][t]));
                maxChange = max(maxChange, fabs(u[i][t] - old_u[i][t]));
            }
        }

        if (maxChange < eps) {
            cout << "Converged at iter = " << iter << ", maxChange = " << maxChange << endl;
            break;
        }
    } // end iterations

    // Profits calculation (as in C#)
    vector<double> profits(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double profit = pow(data.discount, T) * (data.maxCompanyValues[i] - data.liquidity * c[i][T]);
        for (int t = 0; t < T; ++t) {
            double netProfit = data.maxProductCost - c[i][t];
            for (int j = 0; j < n; ++j) netProfit -= u[j][t];
            netProfit *= u[i][t];

            double investmentExpenses = (data.conversion * y[i][t] * y[i][t]) / 2.0;

            double cooperationCosts = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                cooperationCosts += data.cooperationCosts[i][j] * g[i][j][t] * g[j][i][t];
            }

            profit += (netProfit - investmentExpenses - cooperationCosts) * pow(data.discount, t);
        }
        profits[i] = profit;
    }

    // Final adjustments: u[:,T] = 0 and if g[i,j,t]==0 then g[j,i,t]=0
    for (int i = 0; i < n; ++i) {
        u[i][T] = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int t = 0; t <= T; ++t) {
                if (i == j) continue;
                if (g[i][j][t] == 0) g[j][i][t] = 0;
            }
        }
    }

    // Package output
    OutputData out;
    out.y = y;
    out.u = u;
    out.time = T;
    out.profits = profits;
    out.g = g;
    out.num = n;
    return out;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    InputData data;
    cout << "Enter number of firms (num): ";
    cin >> data.num;
    cout << "Enter time horizon (T): ";
    cin >> data.time;
    cout << "Enter discount (rho): ";
    cin >> data.discount;
    cout << "Enter liquidity (eta): ";
    cin >> data.liquidity;

    data.unitCosts.assign(data.num, vector<double>(data.time + 1, 0.0));
    cout << "Enter initial unit costs c_i(0) for each firm:" << endl;
    for (int i = 0; i < data.num; ++i) {
        double v;
        cout << " c: ";
        cin >> v;
        data.unitCosts[i][0] = v;
    }

    cout << "Enter personalInvestments (alpha): ";
    cin >> data.personalInvestments;
    cout << "Enter conversion (epsilon): ";
    cin >> data.conversion;
    cout << "Enter maxProductCost (p): ";
    cin >> data.maxProductCost;

    data.cooperationCosts.assign(data.num, vector<double>(data.num, 0.0));
    cout << "Enter cooperationCosts matrix (num x num):" << endl;
    for (int i = 0; i < data.num; ++i) {
        for (int j = 0; j < data.num; ++j) {
            cout << " coop[" << i + 1 << "][" << j + 1 << "]: ";
            cin >> data.cooperationCosts[i][j];
        }
    }

    cout << "Enter friendsInvestments (beta): ";
    cin >> data.friendsInvestments;
    cout << "Enter notFriendsInvestments (gamma): ";
    cin >> data.notFriendsInvestments;

    data.maxCompanyValues.assign(data.num, 0.0);
    cout << "Enter maxCompanyValues (eta_i) for each firm:" << endl;
    for (int i = 0; i < data.num; ++i) {
        cout << " maxCompanyValue[" << i + 1 << "]: ";
        cin >> data.maxCompanyValues[i];
    }

    cout << "Enter aging (delta): ";
    cin >> data.aging;

    // Solve
    OutputData out = Solve(data);

    // Print summary
    cout << fixed << setprecision(6);
    cout << "\n=== Results ===\n";
    for (int i = 0; i < out.num; ++i) {
        cout << "\n--- Firm " << i + 1 << " ---\n";
        cout << "u(t): ";
        for (int t = 0; t <= out.time; ++t) cout << out.u[i][t] << (t == out.time ? "\n" : ", ");
        cout << "y(t): ";
        for (int t = 0; t <= out.time; ++t) cout << out.y[i][t] << (t == out.time ? "\n" : ", ");
        cout << "Profit (discounted total): " << out.profits[i] << "\n";
    }

    cout << "\nAdjacency matrices g[t] (i->j):\n";
    for (int t = 0; t <= out.time; ++t) {
        cout << "Time t = " << t << ":\n";
        for (int i = 0; i < out.num; ++i) {
            for (int j = 0; j < out.num; ++j) {
                cout << out.g[i][j][t] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }

    return 0;
}
