
#define EPS	(1/1024)

#define ZMAX (20)

/*
Displacement using a z-buffer texture, generated using a surface shader with:
	output varying float Zdist = 0;
	point Pcam = transform("camera", P);
	Zdist = zcomp(Pcam);
Render it with a screen-filling geometry (eg plane or sphere...)
*/

displacement
z_displ(
	string displ_tex = "";
	string mask_tex = ""; )
{
	if( displ_tex != "" )
	{
		point PNDC = transform("NDC", P);
		float ts = xcomp(PNDC);
		float tt = ycomp(PNDC);
		float zdist = ZMAX * float texture(displ_tex, ts, tt);


		float mask = 1;
		if (mask_tex != "") {
			mask = float texture(mask_tex, xcomp(PNDC), ycomp(PNDC));
		}

		if (mask >= 0.90) {
			P = E + normalize(I) * zdist;

			//TODO: add noise based on temperature (?)

		}

		N = calculatenormal(P);
	}
}

