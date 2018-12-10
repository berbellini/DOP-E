#include <stdio.h>
#include <math.h>


sdr_to_mij(s,d,r,m)
float s,d,r,m[3][3];
{
      float d2r,sd,cd,s2d,c2d,ss,cs,c2s,s2s,sr,cr;

      d2r = acos(-1) / 180;
					
      sd = sin(d2r * d);
      cd = cos(d2r * d);
      s2d = sin(d2r * 2 * d);
      c2d = cos(d2r * 2 * d);
      
      ss = sin(d2r * s);
      cs = cos(d2r * s);
      s2s = sin(d2r * 2 * s);
      c2s = cos(d2r * 2 * s);
      
      sr = sin(d2r * r);
      cr = cos(d2r * r);
      
      m[0][0] = -sd*cr*s2s - s2d*sr*ss*ss;
      m[1][1] =  sd*cr*s2s - s2d*sr*cs*cs;
			m[2][2] =  s2d*sr;
      m[0][1] =  sd*cr*c2s + 0.5*s2d*sr*s2s;
      m[0][2] = -cd*cr*cs  - c2d*sr*ss;
      m[1][2] = -cd*cr*ss  + c2d*sr*cs;
      
			m[1][0] = m[0][1];
			m[2][1] = m[1][2];
			m[2][0] = m[0][2];

      return;
}

ar_to_hrv_mij(m)
float m[3][3];
{
      float tmp[3][3];
      int x=0,y=1,z=2,r=0,t=1,f=2, i,j;
      
      tmp[r][r] =  m[z][z];
			tmp[r][t] =  m[z][x];
			tmp[r][f] = -m[z][y];
			
			tmp[t][r] =  tmp[r][t];
			tmp[t][t] =  m[x][x];
			tmp[t][f] = -m[x][y];
			
			tmp[f][r] =  tmp[r][f];
			tmp[f][t] =  tmp[t][f];
			tmp[f][f] =  m[y][y];
			
			for(i=0;i<3;i++)
			 for(j=0;j<3;j++)
				 m[i][j] = tmp[i][j];
			
      return;
}

hrv_to_ar_mij(m)
float m[3][3];
{
      float tmp[3][3];
      int x=0,y=1,z=2,r=0,t=1,f=2, i,j;
			
      tmp[x][x] =  m[2][2];
			tmp[x][y] = -m[t][f];
			tmp[x][z] =  m[t][r];
			
			tmp[y][x] =  tmp[x][y];
			tmp[y][y] =  m[f][f];
			tmp[y][z] = -m[f][r];
			
			tmp[z][x] =  tmp[x][z];
			tmp[z][y] =  tmp[y][z];
			tmp[z][z] =  m[r][r];
			
			for(i=0;i<3;i++)
			 for(j=0;j<3;j++)
				 m[i][j] = tmp[i][j];
			
      return;
}
