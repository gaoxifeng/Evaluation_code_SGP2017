/*
   Copyright (C)  2000    Daniel A. Atkinson  <DanAtk@aol.com>
   Copyright (C)  2004    Ivano Primi  <ivprimi@libero.it>

   This file is part of the HPA Library.

   The HPA Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The HPA Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the HPA Library; if not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301 USA.
*/

/* Constants file for XDIM=23
*/

const int xItt_div = 3;
const int xK_tanh = 7;
const int xMS_exp = 57;
const int xMS_hyp = 65;
const int xMS_trg = 75;

const struct xpr xPi4 = { {
    0x3ffe, 0xc90f, 0xdaa2, 0x2168,
    0xc234, 0xc4c6, 0x628b, 0x80dc,
    0x1cd1, 0x2902, 0x4e08, 0x8a67,
    0xcc74, 0x020b, 0xbea6, 0x3b13,
    0x9b22, 0x514a, 0x0879, 0x8e34,
    0x04dd, 0xef95, 0x19b3, 0xcd3a
  }
};
const struct xpr xPi2 = { {
    0x3fff, 0xc90f, 0xdaa2, 0x2168,
    0xc234, 0xc4c6, 0x628b, 0x80dc,
    0x1cd1, 0x2902, 0x4e08, 0x8a67,
    0xcc74, 0x020b, 0xbea6, 0x3b13,
    0x9b22, 0x514a, 0x0879, 0x8e34,
    0x04dd, 0xef95, 0x19b3, 0xcd3a
  }
};
const struct xpr xPi = { {
    0x4000, 0xc90f, 0xdaa2, 0x2168,
    0xc234, 0xc4c6, 0x628b, 0x80dc,
    0x1cd1, 0x2902, 0x4e08, 0x8a67,
    0xcc74, 0x020b, 0xbea6, 0x3b13,
    0x9b22, 0x514a, 0x0879, 0x8e34,
    0x04dd, 0xef95, 0x19b3, 0xcd3a
  }
};
const struct xpr xEe = { {
    0x4000, 0xadf8, 0x5458, 0xa2bb,
    0x4a9a, 0xafdc, 0x5620, 0x273d,
    0x3cf1, 0xd8b9, 0xc583, 0xce2d,
    0x3695, 0xa9e1, 0x3641, 0x1464,
    0x33fb, 0xcc93, 0x9dce, 0x249b,
    0x3ef9, 0x7d2f, 0xe363, 0x630c
  }
};
const struct xpr xLn2 = { {
    0x3ffe, 0xb172, 0x17f7, 0xd1cf,
    0x79ab, 0xc9e3, 0xb398, 0x03f2,
    0xf6af, 0x40f3, 0x4326, 0x7298,
    0xb62d, 0x8a0d, 0x175b, 0x8baa,
    0xfa2b, 0xe7b8, 0x7620, 0x6deb,
    0xac98, 0x5595, 0x52fb, 0x4afa
  }
};
const struct xpr xLn10 = { {
    0x4000, 0x935d, 0x8ddd, 0xaaa8,
    0xac16, 0xea56, 0xd62b, 0x82d3,
    0x0a28, 0xe28f, 0xecf9, 0xda5d,
    0xf90e, 0x83c6, 0x1e82, 0x01f0,
    0x2d72, 0x962f, 0x02d7, 0xb1a8,
    0x105c, 0xcc70, 0xcbc0, 0x2c5f
  }
};
const struct xpr xSqrt2 = { {
    0x3fff, 0xb504, 0xf333, 0xf9de,
    0x6484, 0x597d, 0x89b3, 0x754a,
    0xbe9f, 0x1d6f, 0x60ba, 0x893b,
    0xa84c, 0xed17, 0xac85, 0x8333,
    0x9915, 0x4afc, 0x8304, 0x3ab8,
    0xa2c3, 0xa8b1, 0xfe6f, 0xdc84
  }
};

const struct xpr xLog2_e = { {
    0x3fff, 0xb8aa, 0x3b29, 0x5c17,
    0xf0bb, 0xbe87, 0xfed0, 0x691d,
    0x3e88, 0xeb57, 0x7aa8, 0xdd69,
    0x5a58, 0x8b25, 0x166c, 0xd1a1,
    0x3247, 0xde1c, 0x43f7, 0x5517,
    0x6cd6, 0x24d9, 0x2f75, 0xc16c
  }
};
const struct xpr xLog2_10 = { {
    0x4000, 0xd49a, 0x784b, 0xcd1b,
    0x8afe, 0x492b, 0xf6ff, 0x4daf,
    0xdb4c, 0xd96c, 0x55fe, 0x37b3,
    0xad4e, 0x91b6, 0xac80, 0x82e7,
    0x859d, 0x0665, 0x0fde, 0x9dd5,
    0x1f3a, 0x3e24, 0xbeab, 0x63ad
  }
};
const struct xpr xLog10_e = { {
    0x3ffd, 0xde5b, 0xd8a9, 0x3728,
    0x7195, 0x355b, 0xaaaf, 0xad33,
    0xdc32, 0x3ee3, 0x4602, 0x45c9,
    0xa202, 0x3a3f, 0x2d44, 0xf78e,
    0xa53c, 0x7542, 0x4efa, 0x1402,
    0xf3f2, 0x9223, 0x5592, 0xc646
  }
};

const struct xpr xRndcorr = { {
    0x3ffe, 0x8000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x03ee
  }
};

const struct xpr xFixcorr = { {
    0x3e97, 0xc000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000,
    0x0000, 0x0000, 0x0000, 0x0000
  }
};

const struct xpr xNaN = { {
    0x0000, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff
  }
};
