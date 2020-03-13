double I1(double v1, double v2, double a, double R0){
        double dv=v1-v2;
        return dv/(R0+a*fabs(dv));
}
double partial_v1(double v1, double v2, double a, double R0){
        double dv=v1-v2, x=(R0+a*fabs(dv));
        return (x-a*dv)/(x*x);
}
double partial_v2(double v1, double v2, double a, double R0){
        double dv=v1-v2, x=(R0+a*fabs(dv));
        return (-x+a*dv)/(x*x);
}
void updateF1 (VEC &F, VEC v, double V, double a, double R0){
        F[0] = I1(v[0], v[1], a, R0) + I1(v[0], v[3], a, R0) + I1(v[0],V, a, R0);
        F[1] = I1(v[1], v[0], a, R0) + I1(v[1], 0, a, R0);
        F[2] = I1(v[2], v[3], a, R0) + I1(v[2], v[4], a, R0) + I1(v[2], V, a, R0);
        F[3] = I1(v[3], v[0], a, R0) + I1(v[3], v[2], a, R0) + I1(v[3], v[5], a, R0) + I1(v[3], 0, a, R0);
        F[4] = I1(v[4], v[2], a, R0) + I1(v[4], v[5], a, R0);
        F[5] = I1(v[5], v[3], a, R0) + I1(v[5], v[4], a, R0) + I1(v[5], v[6], a, R0);
        F[6] = I1(v[6], v[5], a, R0) + I1(v[6], 0, a, R0);
        return;
}
// Q2
double I2(double v1, double v2, double T, double k, double R0){
        double dv=v1-v2;
        return dv/(R0+k*T);
}
double temp(double v1, double v2, double T, double b, double k, double R0){
        double dv=v1-v2;
        return T - b*dv*dv/(R0+k*T);
}
void stamping1(MAT &J, VEC v, int i, int j, int t, double b, double k, double R0){
        double dv=v[i]-v[j], Rt=(R0+k*v[t]);
        J[i][i] += 1/Rt;
        J[j][j] += 1/Rt;
        J[i][j] -= 1/Rt;
        J[j][i] -= 1/Rt;
        J[i][t] = -(dv/(Rt*Rt))*k;
        J[j][t] = (dv/(Rt*Rt))*k;
        J[t][i] = -2*b*dv/Rt;
        J[t][j] = 2*b*dv/Rt;
        J[t][t] = 1 + (b*dv*dv/(Rt*Rt))*k;^M

        return;
}^M
void stamping2(MAT &J, VEC v, int i, double vj, int t, double b, double k, double R0){^M
        double dv=v[i]-vj, Rt=(R0+k*v[t]);^M
        J[i][i] += 1/Rt;^M
        J[i][t] = -(dv/(Rt*Rt))*k;^M
        J[t][i] = -2*b*dv/Rt;^M
        J[t][t] = 1 + (b*dv*dv/(Rt*Rt))*k;^M
        ^M
        return;^M
}
void updateF2 (VEC &F, VEC v, double V, double b, double k, double R0){
        F[0] = I2(v[0], v[1], v[8], k, R0) + I2(v[0], v[3], v[10], k, R0) + I2(v[0],V, v[7], k, R0);
        F[1] = I2(v[1], v[0], v[8], k, R0) + I2(v[1], 0, v[11], k, R0);
        F[2] = I2(v[2], v[3], v[12], k, R0) + I2(v[2], v[4], v[14], k, R0) + I2(v[2], V, v[9], k, R0);
        F[3] = I2(v[3], v[0], v[10], k, R0) + I2(v[3], v[2], v[12], k, R0) + I2(v[3], v[5], v[15], k, R0) + I2(v[3], 0, v[13], k, R0);
        F[4] = I2(v[4], v[2], v[14], k, R0) + I2(v[4], v[5], v[17], k, R0);
        F[5] = I2(v[5], v[3], v[15], k, R0) + I2(v[5], v[4], v[17], k, R0) + I2(v[5], v[6], v[18], k, R0);
        F[6] = I2(v[6], v[5], v[18], k, R0) + I2(v[6], 0, v[16], k, R0);
        F[7] = temp(v[0], V, v[7], b, k, R0);
        F[8] = temp(v[0], v[1], v[8], b, k, R0);
        F[9] = temp(v[2], V, v[9], b, k, R0);
        F[10] = temp(v[0], v[3], v[10], b, k, R0);
        F[11] = temp(v[1], 0, v[11], b, k, R0);
        F[12] = temp(v[2], v[3], v[12], b, k, R0);
        F[13] = temp(v[3], 0, v[13], b, k, R0);
        F[14] = temp(v[2], v[4], v[14], b, k, R0);
        F[15] = temp(v[3], v[5], v[15], b, k, R0);
        F[16] = temp(v[6], 0, v[16], b, k, R0);
        F[17] = temp(v[4], v[5], v[17], b, k, R0);
        F[18] = temp(v[5], v[6], v[18], b, k, R0);
        return;
}^M
int main()^M
{^M
        int i, count, n=19; // n=7 to solve Q1, n=19 to solve Q2        
        double a=0.1, R0=1, V=0, err=0.000000001, b=1, k=1, T8, T13, T18;
        MAT J(n), LU(n);
        VEC v(n), dv(n), F(n), y(n);

        freopen("o2.txt", "w", stdout);
        // initialize
        for (i=0; i<n; i++) v[i] = 0;
        // Q1
/*
        updateF1 (F, v, V, a, R0);
        for (;V<=5; V+=0.1) {
                count = 0;
                err = 1 + err;
                while (err >0.000000001) {
                        // Jacobian
                        J[0][0] = partial_v1(v[0], v[1], a, R0) + partial_v1(v[0], v[3], a, R0) + partial_v1(v[0], V, a, R0);
                        J[0][1] = partial_v2(v[0], v[1], a, R0);
                        J[0][3] = partial_v2(v[0], v[3], a, R0);
                        J[1][1] = partial_v1(v[1], v[0], a, R0) + partial_v1(v[1], 0, a, R0);
                        J[1][0] = partial_v2(v[1], v[0], a, R0);
                        J[2][2] = partial_v1(v[2], v[3], a, R0) + partial_v1(v[2], v[4], a, R0) + partial_v1(v[2], V, a, R0);
                        J[2][3] = partial_v2(v[2], v[3], a, R0);
                        J[2][4] = partial_v2(v[2], v[4], a, R0);
                        J[3][3] = partial_v1(v[3], v[0], a, R0) + partial_v1(v[3], v[2], a, R0) + partial_v1(v[3], v[5], a, R0) + partial_v1(v[3], 0, a, R0);
                        J[3][0] = partial_v2(v[3], v[0], a, R0);
                        J[3][2] = partial_v2(v[3], v[2], a, R0);
                        J[3][5] = partial_v2(v[3], v[5], a, R0);
                        J[4][4] = partial_v1(v[4], v[2], a, R0) + partial_v1(v[4], v[5], a, R0);
                        J[4][2] = partial_v2(v[4], v[2], a, R0);
                        J[4][5] = partial_v2(v[4], v[5], a, R0);        
                        J[5][5] = partial_v1(v[5], v[3], a, R0) + partial_v1(v[5], v[4], a, R0) + partial_v1(v[5], v[6], a, R0);
                        J[5][3] = partial_v2(v[5], v[3], a, R0);
                        J[5][4] = partial_v2(v[5], v[4], a, R0);
                        J[5][6] = partial_v2(v[5], v[6], a, R0);
                        J[6][6] = partial_v1(v[6], v[5], a, R0) + partial_v1(v[6], 0, a, R0);
                        J[6][5] = partial_v2(v[6], v[5], a, R0);
                        // LU decomposition
                        LU = J;
                         LU = J;
                        LU = luFact(LU);^M
                        y = fwdSubs(LU, -F);^M
                        dv = bckSubs(LU, y);
                        v += dv;
                        // update F
                        updateF1 (F, v, V, a, R0);
                        err = norm1(F);
                        count++;
                }
                cout << V << " " << count;
                cout << " I " << I1(v[1], 0, a, R0) + I1(v[3], 0, a, R0) + I1(v[6], 0, a, R0);
                cout << " Ir2 " << I1(v[0], v[1], a, R0);
                cout << " Ir7 " << I1(v[3], 0, a, R0);
                cout << " Ir12 " << I1(v[5], v[6], a, R0) << endl;
        }
*/
        // Q2
/**/
        updateF2 (F, v, V, b, k, R0);^M
        T8 = 0; T13 = 0; T18 = 0;
        for (;V<=5; V+=0.1) {
                count = 0;
                err = 1 + err;
                while (err >0.000000001) {
                        // initialization^M
                        for (i=0; i<7; i++) J[i][i]=0;
                        // Jacobian^M
                        J[0][1] = 0;    J[1][0] = 0;^M
                        J[0][3] = 0;    J[3][0] = 0;^M
                        J[2][3] = 0;    J[3][2] = 0;^M
                        J[2][4] = 0;    J[4][2] = 0;^M
                        J[3][5] = 0;    J[5][3] = 0;^M
                        J[4][5] = 0;    J[5][4] = 0;^M
                        J[5][6] = 0;    J[6][5] = 0;
                        stamping2(J, v, 0, V, 7, b, k, R0);^M
                        stamping1(J, v, 0, 1, 8, b, k, R0);
                        stamping2(J, v, 2, V, 9, b, k, R0);^M
                        stamping1(J, v, 0, 3, 10, b, k, R0);^M
                        stamping2(J, v, 1, 0, 11, b, k, R0);^M
                        stamping1(J, v, 2, 3, 12, b, k, R0);^M
                        stamping2(J, v, 3, 0, 13, b, k, R0);^M
                        stamping1(J, v, 2, 4, 14, b, k, R0);^M
                        stamping1(J, v, 3, 5, 15, b, k, R0);^M
                        stamping2(J, v, 6, 0, 16, b, k, R0);^M
                        stamping1(J, v, 4, 5, 17, b, k, R0);^M
                        stamping1(J, v, 5, 6, 18, b, k, R0);
                        // LU decomposition
                        LU = J;
                        LU = luFact(LU);^M
                        y = fwdSubs(LU, -F);^M
                        dv = bckSubs(LU, y);
                        v += dv;
                        // update F
                        updateF2 (F, v, V, b, k, R0);
                        //F.printv();
                        err = norm1(F);
                        count++;}
                cout << V << " " << count;^M
                cout << " I " << I2(v[1], 0, v[11], k, R0) + I2(v[3], 0, v[13], k, R0) + I2(v[6], 0, v[16], k, R0);^M
                cout << " Tr2 " << v[8]-T8 << " Tr7 " << v[13]-T13 << " Tr12 " << v[18]-T18 << endl;^M
                T8 = v[8]; T13 = v[13]; T18 = v[18];
        }
/**/
        return 0;^M
}
