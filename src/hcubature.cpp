#include "gauss.h"
#include <utility>


//multivariate integration routines

int choose(int n, int k) {
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    int result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

/** [combination c n p x]
 * get the [x]th lexicographically ordered set of [p] elements in [n]
 * output is in [c], and should be sizeof(int)*[p]
 * "Algorithm 515: Generation of a Vector from the Lexicographical Index"; Buckles, B. P., and Lybanon, M. ACM Transactions on Mathematical Software, Vol. 3, No. 2, June 1977.
 * User lucaroni from https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
 */

void combination(int* c, int n, int p, int x) {
    int i, r, k = 0;
    for (i = 0; i < p - 1; i++) {
        c[i] = (i != 0) ? c[i - 1] : 0;
        do {
            c[i]++;
            r = choose(n - c[i], p - (i + 1));
            k = k + r;
        } while (k < x);
        k = k - r;
    }
    if (p > 1) c[p - 1] = c[p - 2] + x - k; else c[0] = x;
}

void combos(int k, double lambda, int n, std::vector<std::vector<double>>& p) {
    std::vector<int> c(k);
    int cnk = choose(n, k) + 1;
    for (int i = 1; i != cnk; i++) {
        std::vector<double> temp(n, 0.0);
        combination(c.data(), n, k, i);
        for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
        p.push_back(temp);
    }
}

void increment(std::vector<bool>& index, int k, double lambda, int n, int* c, std::vector<double>& temp) {
    // temp size n, all elements initially zero
    if (index.size() == 0) {
        index.push_back(false);
        for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
        return;
    }
    int first_zero = 0;
    while ((first_zero < static_cast<int>(index.size())) && index[first_zero]) first_zero++;
    if (first_zero == static_cast<int>(index.size())) {
        index.flip();
        for (int j = 0; j != static_cast<int>(index.size()); j++) temp[c[j] - 1] *= -1;
        index.push_back(true);
        temp[c[index.size() - 1] - 1] = -lambda;
    }
    else {
        int fzp1 = first_zero + 1;
        for (int i = 0; i != fzp1; i++) {
            index[i] = !index[i];
            temp[c[i] - 1] *= -1;
        }
    }
}

void signcombos(int k, double lambda, int n, std::vector<std::vector<double>>& p) {
    std::vector<int> c(k);
    int cnk = choose(n, k) + 1;
    for (int i = 1; i != cnk; i++) {
        std::vector<double> temp(n, 0.0);
        combination(c.data(), n, k, i);
        std::vector<bool> index; index.clear();
        int p2k = std::pow(2, k);
        for (int j = 0; j != p2k; j++) {
            increment(index, k, lambda, n, c.data(), temp);
            p.push_back(temp);
        }
    }
}

void gauss_kronrod(double a, double b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p,
    unsigned fdim, double* retval)) {
    double c = 0.5 * (a + b);
    double delta = 0.5 * (b - a);
    double f0;
    integrand(1, &c, pars, 1, &f0);
    double I = f0 * wd7[7], Idash = f0 * gwd7[3];
    for (int i = 0; i != 7; i++) {
        double deltax = delta * xd7[i], cp = c + deltax, cm = c - deltax;
        double  fx;
        integrand(1, &cp, pars, 1, &fx);
        double temp;
        integrand(1, &cm, pars, 1, &temp);
        fx += temp;
        I += fx * wd7[i];
        if (i % 2 == 1) Idash += fx * gwd7[i / 2];
    }
    double V = fabs(delta);
    I *= V;
    Idash *= V;
    out.result = I;
    out.err = fabs(I - Idash);
}

void make_GenzMalik(int n, GenzMalik& g) {
    double l4 = sqrt(9 * 1.0 / 10);
    double l2 = sqrt(9 * 1.0 / 70);
    double l3 = l4;
    double l5 = sqrt(9 * 1.0 / 19);

    int twopn = pow(2, n);

    g.w[0] = twopn * ((12824 - 9120 * n + 400 * n * n) * 1.0 / 19683);
    g.w[1] = twopn * (980.0 / 6561);
    g.w[2] = twopn * ((1820 - 400 * n) * 1.0 / 19683);
    g.w[3] = twopn * (200.0 / 19683);
    g.w[4] = 6859.0 / 19683;
    g.wd[3] = twopn * (25.0 / 729);
    g.wd[2] = twopn * ((265 - 100 * n) * 1.0 / 1458);
    g.wd[1] = twopn * (245.0 / 486);
    g.wd[0] = twopn * ((729 - 950 * n + 50 * n * n) * 1.0 / 729);

    combos(1, l2, n, g.p[0]);
    combos(1, l3, n, g.p[1]);
    signcombos(2, l4, n, g.p[2]);
    signcombos(n, l5, n, g.p[3]);
}

void clean_GenzMalik(GenzMalik& g) {
    for (int j = 0; j != 4; j++) {
        int gpjs = g.p[j].size();
        for (int i = 0; i != gpjs; i++) g.p[j][i].clear();
    }
}

void integrate_GenzMalik(GenzMalik g, int n, const double* a, const double* b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval)) {
    std::vector<double> c(n);
    std::vector<double> deltac(n);

    for (int i = 0; i != n; i++) c[i] = (a[i] + b[i]) / 2;
    for (int i = 0; i != n; i++) deltac[i] = fabs(b[i] - a[i]) / 2;
    double v = 1.0;
    for (int i = 0; i != n; i++) v *= deltac[i];

    if (v == 0.0) {
        out.err = 0.0;
        out.result = 0.0;
        out.kdivide = 0;
        return;
    }

    double f1;
    integrand(n, c.data(), pars, 1, &f1);
    double f2 = 0.0, f3 = 0.0;
    double twelvef1 = 12 * f1;

    double maxdivdiff = 0.0;
    std::vector<double> divdiff(n);
    std::vector<double> p2(n);
    std::vector<double> p3(n);
    std::vector<double> cc(n);

    for (int i = 0; i != n; i++) {

        for (int j = 0; j != n; j++) p2[j] = deltac[j] * g.p[0][i][j];

        for (int j = 0; j != n; j++) cc[j] = c[j] + p2[j];
        double f2i;
        integrand(n, cc.data(), pars, 1, &f2i);
        for (int j = 0; j != n; j++) cc[j] = c[j] - p2[j];
        double temp;
        integrand(n, cc.data(), pars, 1, &temp);
        f2i += temp;


        for (int j = 0; j != n; j++) p3[j] = deltac[j] * g.p[1][i][j];
        for (int j = 0; j != n; j++) cc[j] = c[j] + p3[j];
        double f3i;
        integrand(n, cc.data(), pars, 1, &f3i);
        for (int j = 0; j != n; j++) cc[j] = c[j] - p3[j];
        integrand(n, cc.data(), pars, 1, &temp);
        f3i += temp;
        f2 += f2i;
        f3 += f3i;
        divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);

    }
    std::vector<double> p4(n);
    double f4 = 0.0;
    int gp2s = g.p[2].size(), gp3s = g.p[3].size();
    for (int i = 0; i != gp2s; i++) {

        for (int j = 0; j != n; j++) p4[j] = deltac[j] * g.p[2][i][j];
        for (int j = 0; j != n; j++) cc[j] = c[j] + p4[j];
        double temp;
        integrand(n, cc.data(), pars, 1, &temp);
        f4 += temp;
    }
    double f5 = 0.0;
    std::vector<double> p5(n);
    for (int i = 0; i != gp3s; i++) {

        for (int j = 0; j != n; j++) p5[j] = deltac[j] * g.p[3][i][j];

        for (int j = 0; j != n; j++) cc[j] = c[j] + p5[j];
        double temp;
        integrand(n, cc.data(), pars, 1, &temp);
        f5 += temp;
    }
    double I = v * (g.w[0] * f1 + g.w[1] * f2 + g.w[2] * f3 + g.w[3] * f4 + g.w[4] * f5);
    double Idash = v * (g.wd[0] * f1 + g.wd[1] * f2 + g.wd[2] * f3 + g.wd[3] * f4);
    double E = fabs(I - Idash);

    int kdivide = 0;
    double deltaf = E / (pow(10, n) * v);
    for (int i = 0; i != n; i++) {
        double delta = divdiff[i] - maxdivdiff;
        if (delta > deltaf) {
            kdivide = i;
            maxdivdiff = divdiff[i];
        }
        else if ((fabs(delta) <= deltaf) && (deltac[i] > deltac[kdivide])) kdivide = i;
    }
    out.result = I;
    out.err = E;
    out.kdivide = kdivide;
}

class Box {
public:
	Box(std::vector<double> a, std::vector<double> b, double I, double err, int kdivide) : a(std::move(a)), b(std::move(b)), I(I), E(err), kdiv(kdivide) {};
	bool operator<(const Box& box) const { return E < box.E; }
	std::vector<double> a;
	std::vector<double> b;
	double I;
	double E;
	int kdiv;
};

Box make_box(int n, const double* a, const double* b, one_d out) {
	return Box(std::vector<double>(a, a + n), std::vector<double>(b, b + n),
		out.result, out.err, out.kdivide);
}

int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval), void* pars, unsigned n, const double* a, const double* b,
    size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err) {

    one_d out;
    GenzMalik g;

    if (n == 1) gauss_kronrod(a[0], b[0], out, pars, integrand);
    else {
        make_GenzMalik(n, g);
        integrate_GenzMalik(g, n, a, b, out, pars, integrand);
    }
    int numevals = (n == 1) ? 15 : 1 + 4 * n + 2 * n * (n - 1) + pow(2, n);
    int evals_per_box = numevals;
    // int kdiv = out.kdivide;
    err[0] = out.err;
    val[0] = out.result;
    // convergence test
    if ((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval!=0) && (numevals >= maxEval))) {
//        std::cout << numevals << std::endl;
        return 0;
    }

    std::priority_queue<Box> ms;
    ms.push(make_box(n, a, b, out));

    while (true) {
        Box box = ms.top();
        ms.pop();
        // split along dimension kdiv
        double w = (box.b[box.kdiv] - box.a[box.kdiv]) / 2;
        std::vector<double> ma = box.a;
        ma[box.kdiv] += w;
        std::vector<double> mb = box.b;
        mb[box.kdiv] -= w;

        if (n == 1) gauss_kronrod(ma[0], box.b[0], out, pars, integrand);
        else {
            integrate_GenzMalik(g, n, ma.data(), box.b.data(), out, pars, integrand);
        }
        Box box1 = make_box(n, ma.data(), box.b.data(), out);
        ms.push(box1);

        if (n == 1) gauss_kronrod(box.a[0], mb[0], out, pars, integrand);
        else {
            integrate_GenzMalik(g, n, box.a.data(), mb.data(), out, pars, integrand);
        }
        Box box2 = make_box(n, box.a.data(), mb.data(), out);
        ms.push(box2);
        val[0] += box1.I + box2.I - box.I;
        err[0] += box1.E + box2.E - box.E;
        numevals += 2 * evals_per_box;
        if (((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval != 0) && (numevals >= maxEval))) || !(std::isfinite(val[0])) ) {
            break;
        }
    }
    val[0] = 0.0;
    err[0] = 0.0;

    while (!ms.empty()) {
        Box box = ms.top();
        val[0] += box.I;
        err[0] += box.E;
        ms.pop();
    }
    clean_GenzMalik(g);
    return 0;
}
