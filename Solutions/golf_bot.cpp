#include <bits/stdc++.h>
using namespace std;
const double PI{acos(-1.0)}, EPS{1e-6};

void fft(vector<complex<double>> &xs, bool invert = false){
    int N = (int)xs.size();
    if (N == 1) return;

    vector<complex<double>> es(N / 2), os(N / 2);
    for (int i = 0; i < N / 2; ++i) es[i] = xs[2 * i];
    for (int i = 0; i < N / 2; ++i) os[i] = xs[2 * i + 1];

    fft(es, invert);
    fft(os, invert);

    auto signal = (invert ? 1 : -1);
    auto theta = 2 * signal * PI / N;
    complex<double> S{1}, S1{cos(theta), sin(theta)};
    for (int i = 0; i < N / 2; ++i){
        xs[i] = (es[i] + S * os[i]);
        xs[i] /= (invert ? 2 : 1);

        xs[i + N / 2] = (es[i] - S * os[i]);
        xs[i + N / 2] /= (invert ? 2 : 1);

        S *= S1;
    }
}
int main(){
    int N;
    while (cin >> N){
        int ans = 0, M, size = 1;
        vector<int> ki(N);
        for (int i = 0; i < N; ++i) cin >> ki[i];

        cin >> M;

        vector<int> dj(M);

        for (int j = 0; j < M; ++j) cin >> dj[j];
        auto K = *max_element(ki.begin(), ki.end());

        while (size < K + 1) size *= 2;
        size *= 2; // Encontrar um tamanho adequado para os vetores dps da multiplicacao

        vector<complex<double>> xs(size, 0);
        for (auto k : ki) xs[k] = 1; // Inicializa o vetor xs com os coeficientes do primeiro

        xs[0] = 1;
        fft(xs); // Aplica a FFT no vetor

        for (int i = 0; i < size; ++i) xs[i] *= xs[i]; // Calcula o quadrado dos coeficientes transformados

        fft(xs, true); // Aplica a FFT inversa para obter o resultado final

        for (auto d : dj) ans += (d < size and fabs(xs[d].real()) > EPS ? 1 : 0); // Conta o numero de coeficientes nao nulos nas posicoes

        cout << ans << '\n';
    }

    return 0;
}
