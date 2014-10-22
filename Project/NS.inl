updateBoundaryGhosts()
{
	//u_fluxes
	DMDAVecGetArray(uda, qxLocal, &qx);
	DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL);
	DMDAGetInfo();

	//x-faces
	
