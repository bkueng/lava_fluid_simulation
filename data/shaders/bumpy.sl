/******************************************************************************/
/*                                                                            */
/*    Copyright (c)The 3Delight Team.                                         */
/*    All Rights Reserved.                                                    */
/*                                                                            */
/******************************************************************************/

displacement
bumpy(
	float amplitude = 1;
	string texturename = ""; )
{
	if( texturename != "" )
	{
		float amp = amplitude * float texture( texturename, s, t );

		P += amp * normalize( N );
		N = calculatenormal( P );
	}
}
