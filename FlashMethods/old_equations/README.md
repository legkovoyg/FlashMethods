# FlashMethods

Тут лежат сами скрипты, в том виде в котором я их просто сделал тыщу лет назад
Тут нету никаких dto, проверок, вызовов, просто методы которые вызываются в другом месте через:

import PR_equation as PR
import PR_equation_with_peneloux as PR_peneloux
import SRK_equation as SRK
import SRK_with_peneloux_equation as SRK_peneloux

BIPS = None
pres = 25 #МПа
temp = 300 #K
v = PR.PR_Flash(Fluid, BIPs = BIPS)
( W,
Z_v,
Z_l,
x_i,
y_i,
Stable,
m,
enthalpy,
enthalpy_w,
enthalpy_l,
Cp,
Cp_w,
Cp_l,
Cv,
Cv_w,
Cv_l,
volume,
VolumeMy_y,
VolumeMy_x,
density,
density_y,
density_x) = v.vle(pres, temp)
