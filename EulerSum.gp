\\ Alternating multiple zeta values を計算するための pari/gp　プログラム

\\ 反復積分 I(0;vs[1], … , vs[len(vs)]; 1) の計算。成分は 0 もしくは 絶対値が1より大きい数を想定。
\\ \sum_{n=0}^{\infty} a(d,n) z^n := I(0;vs[1], … , vs[d-1]; z) \times 1/(z-vs[d]) としたとき
\\ a(d,n)は次の漸化式を満たす。 x=vs[d]とおく。
\\ n = 0の場合：
\\ a(d,0)は、d=1のとき -1/vs[1]、d>1 のとき x==0 であれば a(d-1)、そうでなければ 0
\\ n>0の場合：
\\ a(1,n) = -vs[1]^(-n-1)
\\ x =0 のとき、a(d,n) = a(d-1,n)/(n+1)
\\ x!=0 のとき、a(d,n) = ( a(d,n-1) - a(d-1,n-1)/n ) / x
\\
\\ \sum_{n=0}^{\infty} a(len(vs), n) / (n+1) を計算することで反復積分が計算される。
iteratedintegral_2(vs)={
local(ret, a, new_a, n, x, ret2, flag);
if(length(vs)==0, return(1));
ret=0;
a = List();
listput(a, -1.0/vs[1]);
for(d=2, length(vs), listput(a, if(vs[d]==0, a[d-1],0.0 )));
n=0;
flag = 0;
while(1,
	ret2 = ret + a[length(vs)] / (n+1);
	\\ 計算精度を保障するために必要な計算回数は微妙な問題だが、とりあえず値の更新が5連続でなければ計算を打ち切ることにする。
	if(n>length(vs) && ret==ret2, flag += 1, flag=0);
	if(flag>=5, break);
	ret = ret2;

	n += 1;
	new_a = a;
	new_a[1] = -1.0*vs[1]^(-n-1);
	for(d=2, length(vs),
		x = vs[d];
		new_a[d] = if(x==0, new_a[d-1]/(n+1), (a[d]-a[d-1]/n)/x  );
	);
	a = new_a;
);
ret
};




\\ 反復積分 I(0;vs[1], … , vs[len(vs)]; 1) の計算。成分は 0 or 1 or -1 を想定
\\ path composition formulaを用いて計算する。
iteratedintegral_1(vs)={
local(ret=0, vs1, vs2 );
for(i=0,length(vs),
	vs1 = List();
	for(j=1, i, listput(vs1, 2*vs[j]));
	vs2 = List();
	forstep(j=length(vs), i+1, -1, listput(vs2, 2*(1-vs[j])));
	ret += (-1)^(length(vs)-i) * iteratedintegral_2( vs1 ) * iteratedintegral_2( vs2 );
);
ret
};



\\ 交代多重ゼータ値の計算。 例えば AMZV([1,2]) = zeta(3), AMZV([-1]) = -log(2)。
\\ pari/gpの組み込み関数zetamultとはindexの順序が逆であることに注意。
AMZV(ks)={
local(seq=List(), eps=1);
for(i=1, length(ks),
	listput(seq,eps);
	for(i=1, abs(ks[i])-1, listput(seq, 0));
	eps *= sign(ks[i]);
);
if(eps==-1,
	for(i=1, length(seq), seq[i] *= -1);
);
(-1)^length(ks) * iteratedintegral_1(seq)
}




test()={
	\\ 全てのindexが正の場合は多重ゼータ値と一致するはずである。
	print(zetamult([4,3,5]) - AMZV([5,3,4]));

	\\ distribution relationの確認。
	print( AMZV([5,3])+AMZV([5,-3])+AMZV([-5,3])+AMZV([-5,-3]) - AMZV([5,3])/64 );

}
test()
