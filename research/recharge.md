# A second climb from a Mersenne seed

**For every integer $m\ge2$, the orbit from $2^m-1$ cannot immediately repeat its first run of $m$ odd shortcut steps after reaching its next odd state.** The result concerns those two consecutive climbs. The orbit can have later climbs of any length permitted by the map.

There is a useful companion construction. Allow a positive odd multiplier on the power of two in the seed, and the second climb can be exponentially longer than the first. For example, **77,671 takes three consecutive odd shortcut steps and one halving to reach 131,071, which then takes seventeen consecutive odd shortcut steps.**

## Precisely what counts as a climb

Let $T$ halve even numbers and send an odd $x$ to $(3x+1)/2$. An odd shortcut step represents two standard steps. Write $v_2(z)$ for the exponent of 2 in a positive integer $z$. If $x=2^m u-1$ with $u$ odd and positive, then for $0\le j\le m$,

$$T^j(x)=3^j2^{m-j}u-1.$$

The inputs at steps $0$ through $m-1$ are odd, and the endpoint is even. Divide out its remaining powers of 2 to reach the next odd state. In particular, with $u=1$ this gives

$$t=v_2(3^m-1),\qquad y=(3^m-1)/2^t.$$

The next climb has exactly $r=v_2(y+1)$ odd shortcut steps by the same identity. Consequently $r\ge m$ would require

$$v_2(3^m+2^t-1)\ge m+t.$$

For odd $m$, $3^m\equiv3\pmod8$, so $t=1$, $y\equiv1\pmod4$, and $r=1$. The seed at $m=1$ is the trivial orbit through 1. At $m=2$, $y=1$ and $r=1<2$. It remains to treat even $m\ge4$.

## The finite certificate for $m<2^{20}$

Write $m=2^a b$ with $b$ odd. Since $3^2=9$ and each subsequent squaring increases $v_2(3^{2^a}-1)$ by one, $t=a+2$. Raising to an odd power preserves this valuation. The same calculation gives, for $b_1>b_2\ge0$,

$$v_2(3^{2^a b_1}-3^{2^a b_2})=t+v_2(b_1-b_2).$$

For fixed $a$ and depth $h$, the congruence $3^{2^a b}\equiv1-2^t\pmod{2^{t+h}}$ therefore has at most one residue class of solutions $b\bmod 2^h$. If $m\ge h$, full recharge would imply this weaker congruence.

For each $a=1,\ldots,19$, the certificate supplies a depth $h$ and a least nonnegative residue $b_0$. The checker verifies the modular power directly and verifies $\lfloor(2^{20}-1)/2^a\rfloor<b_0<2^h$. This excludes every admissible $b$ for that row when $m\ge h$. For example, $a=1$ gives $h=22$ and $b_0=2{,}490{,}967$, whereas the largest possible $b$ is 524,287. The greatest depth across the 19 rows is 22. Direct integer arithmetic settles $3\le m<22$. The verifier needs no trust in how the lifting algorithm discovered its rows. [Machine-readable certificate](../results/recharge.json) · [Generator and checker](../recharge.py)

## A published valuation estimate covers every larger exponent

Use natural logarithms in this section. For even $m\ge2^{20}$ set $d=\lfloor\log_2m\rfloor$. Apply Bennett, Bugeaud and Mignotte's Theorem 6 to the positive rational bases $\alpha_1=9$ and $\alpha_2=(2^t-1)^2$, with exponents $m$ and 1. They are multiplicatively independent: dependence would force $2^t-1$ to be a power of 3, but it is 7 modulo 8. The theorem's common residue order is $g=1$. We may take $E=3$ because $v_2(9-1)=3$ and $v_2(\alpha_2-1)=t+1\ge4$. [*Perfect powers with few binary digits*, Theorem 6, printed p. 528](https://personal.math.ubc.ca/~bennett/BeBuMi-Camb-2012.pdf)

Let $H_1=\ln 9$ and $H_2=2t\ln 2$. They meet both the theorem's logarithmic height and $E\ln 2$ lower bounds. Specializing its estimate to $\Lambda=9^m-(2^t-1)^2$ gives

$$v_2(\Lambda)\le \frac{36.1}{27(\ln2)^4}\,H_1H_2\,
\max\{\ln b'+\ln(3\ln2)+0.4,\ 18\ln2,\ 5\}^2,\qquad
b'=m/H_2+1/H_1.$$

Here $2/3<\ln2<1$, so the leading coefficient is below 7. Also $H_1<3$, $H_2<2t$, and $b'<m$. Since $\ln m<(d+1)\ln2$ and $\ln(3\ln2)+0.4<2$, the maximum is less than $d+3$. Finally $t=v_2(m)+2\le d+2$. Combining these bounds gives

$$v_2(\Lambda)<42(d+3)^3<2^d\le m.$$

The middle inequality starts at $42\cdot23^3=511{,}014<2^{20}$. Each increment of $d$ multiplies the cubic by at most $(24/23)^3<2$, while $2^d$ doubles.

On factoring $\Lambda$, the factor $3^m-(2^t-1)$ is 2 modulo 8. Therefore $v_2(\Lambda)=1+v_2(3^m+2^t-1)$. Full recharge would force this valuation to be at least $1+m+t$, contradicting the upper bound. Together with the finite certificate, this proves the claim. It also bounds the next climb by $O((\log m)^3)$ throughout the even tail.

## How to force a longer second climb

For any $m\ge2$, set $L=2\cdot3^{m-1}$ and $u=(2^L-1)/3^m$. This is a positive odd integer. Indeed $4\equiv1\pmod3$, and cubing a number congruent to 1 modulo 3 introduces another factor of 3 into its difference from 1. Apply this $m-1$ times.

The climb identity at the top of the page now sends $2^m u-1$ to $3^m u-1=2^L-2$ in $m$ odd shortcut steps. One further halving gives the Mersenne seed $2^{L-1}-1$. Its next climb lasts exactly $L-1$ odd shortcut steps. At $m=3$ the first four shortcut states after 77,671 are 116,507, 174,761, 262,142 and 131,071.

Even with $m$ fixed, any positive multiple of $2\cdot3^{m-1}$ works as $L$, giving arbitrarily long second climbs by changing the cofactor $u$. After the second climb and its remaining halvings, the odd state is the odd part of $3^{L-1}-1$. Since $L-1$ is odd, the following climb has exactly one odd shortcut step. Thus the construction gives two arbitrarily long climbs along an actual orbit and identifies exactly where that pattern ends.

## Reproduce

```powershell
python recharge.py --verify results/recharge.json
python -m unittest -v test_recharge
```

Regenerate the finite certificate with `python recharge.py --bits 20 --output results/recharge.json`. The test suite also checks each exponent in the finite range independently, updating $3^m$ modulo $2^{43}$ and checking that the retained bits determine each valuation. The infinite-tail argument uses the cited theorem, independently of the enumeration. This project has not established literature priority for its application here.
