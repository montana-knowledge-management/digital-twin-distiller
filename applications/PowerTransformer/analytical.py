from math import pi
mu_0 = 4*pi*1e-7
def short_circuit_impedance(b_pow, p_num, freq, alpha, turn_v, h, s,
                            r_in, t_in, r_ou, t_ou, g):
    """
    Short-circuit impedance calculation
    b_pow - built-in power  [kVA]
    p_num - phase number [#]
    freq  - frequency  [Hz]
    alpha - ratio of the outer and inner winding
    ff_c  - core filling factor
    turn_v- turn voltage [V]
    h     - height of the inner window [mm]
    s     - width of working window [mm]
    g     . main insulation width [mm]
    """
    p_pow = b_pow / p_num
    imp_con = 4. * pi ** 2. * mu_0 * freq * p_pow / turn_v ** 2. / (h * (1 + alpha) / 2. + 0.32 * s)
    a = r_in * t_in / 3.
    b = r_ou * t_ou / 3.
    c = (r_in + t_in / 2. + g / 2.) * g
    return imp_con * (a + b + c)

b_pow = 6300.
p_num = 3.
freq = 50.
alpha = 0.95
t_in = 42
r_in = 437. / 2.
t_ou = 41
r_ou = 578./2
g = 27.
h = 979.
s = 152.
turn_v = 31.0
sci = round(short_circuit_impedance(b_pow, p_num, freq, alpha, turn_v, h, s, r_in, t_in, r_ou, t_ou, g), 3)
print(sci*100)
