#include <bits/stdc++.h>

using namespace std;

const double PI{ acos(-1.0) };

int reversed(int x, int bits){ //reverte numero de bits em x
    int res = 0;

    for (int i = 0; i < bits; ++i){
        res <<= 1;
        res |= (x & 1);
        x >>= 1;
    }

    return res;
}

void fft(vector<complex<double>>& xs, bool invert = false){
    int N = (int) xs.size();
    if (N == 1) return;

    int bits = 1;

    while ((1 << bits) != N) ++bits; // encontrar o n de bits para representar N

    for (int i = 0; i < N; ++i) {
        auto j = reversed(i, bits);

        if (i < j) swap(xs[i], xs[j]);
    }

    for (int size = 2; size <= N; size *= 2){  // FFT para radices-2 em tempo de decimacao
        auto signal = (invert ? 1 : -1);
        auto theta = 2 * signal * PI / size;
        complex<double> S1{ cos(theta), sin(theta) };

        for (int i = 0; i < N; i += size){
            complex<double> S{ 1 }, k{ invert ? 2.0 : 1.0 };

            for (int j = 0; j < size / 2; ++j){
                auto a{ xs[i + j] }, b{ xs[i + j + size/2] * S };
                xs[i + j] = (a + b) / k;
                xs[i + j + size/2] = (a - b) / k;
                S *= S1;
            }
        }
    }
}

int main(){
    string s;
    cin >> s;
    const string cs{ "abc" };
    int size = 1, m = 0, N = (int) s.size();

    vector<int> is, ms(N + 1, 0);
    // Ajusta o tamanho para ser uma potencia de 2
    while (size < N + 1) size *= 2;
    size *= 2;

    for (auto c : cs){
        vector<complex<double>> xs(size), ys(size);

        for (int i = 0; i < N; ++i) xs[i] = (s[i] == c ? 1 : 0);
        for (int i = 0; i < N; ++i) ys[i] = xs[N - 1 - i];

        fft(xs);
        fft(ys);

        // multiplica elemento a elemento as sequencias xs e ys
        for (int i = 0; i < size; ++i) xs[i] *= ys[i];
        // soma
        fft(xs, true);
        for (int i = 1; i < N; ++i) ms[i] += (int) round(xs[N - 1 + i].real());
    }
    // encontra a posicoes com a maior soma
    for (int i = 1; i <= N; ++i){
        if (ms[i] > m){
            m = ms[i];
            is = vector<int>{ i };
        } else if (ms[i] == m)
            is.push_back(i);
    }

    auto ans = make_pair(m, is);

    cout << ans.first << '\n';

    for (size_t i = 0; i < ans.second.size(); ++i)
        cout << ans.second[i] << (i + 1 == ans.second.size() ? '\n' : ' ');

    return 0;
}