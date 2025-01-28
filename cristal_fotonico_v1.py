import meep as mp
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot_epsilon(sim, filename="geometry.png"):
    """Plot the geometry of the photonic crystal."""
    eps_data = sim.get_array(center=mp.Vector3(), size=sim.cell_size, component=mp.Dielectric)
    plt.imshow(eps_data.T, interpolation='spline36', cmap='binary', origin='lower')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar(label="Epsilon")
    plt.title("Photonic Crystal Geometry")
    plt.savefig(filename)
    plt.close()

def main(args):
    resolution = 20 #resoluçao da simulaçao (pixels/um)

    eps = 13 #constante dielétrica da guia de onda
    w = 1.2 #largura da guia de onda
    r = 0.36 #raio dos buracos
    d = 1.4 #espaçamento do defeito
    N = args.N #quantidade de buracos em cada lado do defeito

    sy = args.sy #tamanho da celula na direçao y
    pad = 2 #espaçamento entre a borda PML e o ultimo buraco
    dpml = 1 #largura do PML

    sx = 2*(pad+dpml+N)+d-1 #tamanho da celula na direçao x levando em conta os elementos da simulaçao

    cell = mp.Vector3(sx,sy,0) #celula computacional

    blk = mp.Block(size=mp.Vector3(mp.inf,w,mp.inf), material=mp.Medium(epsilon=eps)) #guia de onda
    geometry = [blk]

    for i in range(N):
        geometry.append(mp.Cylinder(r, center=mp.Vector3(d/2+i))) #buracos a direita do defeito
        geometry.append(mp.Cylinder(r, center=mp.Vector3(-(d / 2 + i)))) #buracos a esquerda do defeito (simétrico)

    pml_layers = [mp.PML(1.0)] #camada PML

    fcen = args.fcen #frequencia central da fonte (??)
    df = args.df #intervalo de frequencia (??)

    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df), #qual a natureza dessa fonte ?
                     component=mp.Ey, center=mp.Vector3(-0.5*sx+dpml), size=mp.Vector3(0,w))] #posiçao da fonte

    sym = [mp.Mirror(mp.Y, phase=-1)] #nao entendi porque é simétrico no eixo Z e não no eixo Y

    sim = mp.Simulation(cell_size=cell, #caracteristicas da simulaçao
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        sources=src,
                        symmetries=sym,
                        resolution=resolution)

    freg = mp.FluxRegion(center=mp.Vector3(0.5*sx-dpml-0.5), size=mp.Vector3(0,2*w)) #regiao medida do fluxo (??)

    nfreq = 500 #quantidade de frequencias (??)

    trans = sim.add_flux(fcen, df, nfreq, freg) #adicionando o calculo do fluxo a simulaçao

    vol = mp.Volume(mp.Vector3(0), size=mp.Vector3(sx))

    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.during_sources(mp.in_volume(vol, mp.to_appended("hz-slice", mp.at_every(0.4, mp.output_hfield_z)))),
            until_after_sources = mp.stop_when_fields_decayed(50,mp.Ey,mp.Vector3(0.5*sx-dpml-0.5),1e-3))

    sim.display_fluxes(trans)

    plot_epsilon(sim)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type=int, default=3, help='number of holes on either side of defect')
    parser.add_argument('-sy', type=int, default=6, help='size of cell in y direction (perpendicular to wvg.)')
    parser.add_argument('-fcen', type=float, default=0.25, help='pulse center frequency')
    parser.add_argument('-df', type=float, default=0.2, help='pulse frequency width')
    args = parser.parse_args()
    main(args)
