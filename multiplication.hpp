void mult(double *a, double *b, double *c, int n, int m)
{
    int k, l, i, j, t, q, r, s, bl;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    double *pa, *pb, *pc;
    k = n/m;
    l = n - k*m;
    bl = (l!=0 ? k+1 : k);
    for(i=0; i<bl; i++)
        for(j=0; j<bl; j++)
        {//размер Cij v*h
            int v, h;
            v = (i<k ? m : l);
            h = (j<k ? m : l);
            pc = c + i*n*m + j*m;
            for(r=0; r<v; r++)
                for(t=0; t<h; t++)
                    pc[r*n+t] = 0;
            //весь блок мат-цы С попадает в кэш память
            for(s=0; s<bl; s++)
            {
                int ah, v3, h3;
                ah = (s<k ? m : l); //столбцов в A
                pa = a + i*n*m + s*m; //Ais
                pb = b + s*n*m + j*m; //Bsj
                // умножаем Ais Bsj
                v3 = v%3;
                h3 = h%3;
                for(r=0; r<v3; r++)
                {
                    for(t=0; t<h3; t++)
                    {
                        double sum;
                        sum=0;
                        for(q=0; q<ah; q++)
                            sum+=pa[r*n+q]*pb[q*n+t];
                        pc[r*n+t]+=sum;
                    }
                    for(; t<h; t+=3)
                    {
                        s00=0;
                        s01=0;
                        s02=0;
                        for(q=0; q<ah; q++)
                        {
                            s00 += pa[r*n+q]*pb[q*n+t];
                            s01 += pa[r*n+q]*pb[q*n+t+1];
                            s02 += pa[r*n+q]*pb[q*n+t+2];
                        }
                        pc[r*n+t] += s00;
                        pc[r*n+t+1] += s01;
                        pc[r*n+t+2] += s02;
                    }
                }
                for(;r<v; r+=3)
                {
                    for(t=0; t<h3; t++)
                    {
                        s00 = 0;
                        s10 = 0;
                        s20 = 0;
                        for(q=0; q<ah; q++)
                        {
                            s00 += pa[r*n+q]*pb[q*n+t];
                            s10 += pa[(r+1)*n+q]*pb[q*n+t];
                            s20 += pa[(r+2)*n+q]*pb[q*n+t];
                        }
                        pc[r*n+q]+=s00;
                        pc[(r+1)*n+q]+=s10;
                        pc[(r+2)*n+q]+=s20;
                    }
                    for(; t<h; t+=3)
                    {
                        s00=0;
                        s01=0;
                        s02=0;
                        s10=0;
                        s11=0;
                        s12=0;
                        s20=0;
                        s21=0;
                        s22=0;
                        for(q=0; q<ah; q++)
                        {
                            s00 += pa[r*n+q]*pb[q*n+t];
                            s01 += pa[r*n+q]*pb[q*n+t+1];
                            s02 += pa[r*n+q]*pb[q*n+t+2];
                            s10 += pa[(r+1)*n+q]*pb[q*n+t];
                            s11 += pa[(r+1)*n+q]*pb[q*n+t+1];
                            s12 += pa[(r+1)*n+q]*pb[q*n+t+2];
                            s20 += pa[(r+2)*n+q]*pb[q*n+t];
                            s21 += pa[(r+2)*n+q]*pb[q*n+t+1];
                            s22 += pa[(r+2)*n+q]*pb[q*n+t+2];
                        }
                        pc[r*n+t] += s00;
                        pc[r*n+t+1] += s01;
                        pc[r*n+t+2] += s02;
                        pc[(r+1)*n+t] += s10;
                        pc[(r+1)*n+t+1] += s11;
                        pc[(r+1)*n+t+2] += s12;
                        pc[(r+2)*n+t] += s20;
                        pc[(r+2)*n+t+1] += s21;
                        pc[(r+2)*n+t+2] += s22;
                    }
                }
            }
        }
}
