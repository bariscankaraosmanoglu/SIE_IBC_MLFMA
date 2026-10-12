#include "generate_edgeconnection.h"
#include "surface_components.h"

void generate_edgeconnection(){
    generate_surface_components();
    //
    //
    // generate_edgeconnection
    // Generate edge connectivity for the triangles
    int *num_trispernode;
    int size_nodeelements=0;
    int indnode=0;
    int ind0=0;
    int nodes[3]{0,0,0};
    int tri1,tri2;
    int element1node[3]{0,0,0};
    int opengeometry=0;
    double rwg_nodepos1[3]{0,0,0};
    double rwg_nodepos2[3]{0,0,0};
    double temp_node_pos1[3];
    double temp_node_pos2[3];
    double temp_node_pos3[3];
    double temp_vec1[3];
    double temp_vec2[3];
    double temp_crossprod[3];
    double temp_normalize;
    num_trispernode = new int[num_nodes];
    nodetrisindex = new int[num_nodes];

    for(int ind1=0;ind1<(num_nodes);ind1++){
        num_trispernode[ind1]=0;
        nodetrisindex[ind1]=0;
    }
    for(int ind1=0;ind1<(num_tris*3);ind1++){
        num_trispernode[tri_nodes[ind1]]++;
    }
    for(int ind1=0;ind1<num_nodes;ind1++){
        size_nodeelements+=num_trispernode[ind1];
    }
    for(int ind1=1;ind1<num_nodes;ind1++){
        nodetrisindex[ind1]=nodetrisindex[ind1-1]+num_trispernode[ind1-1];
    }
    int *nodetris = new int[size_nodeelements];
    for(int ind1=0;ind1<num_tris;ind1++){
        for(int ind2=0;ind2<3;ind2++){
            indnode=tri_nodes[ind0];
            nodetris[nodetrisindex[indnode]]=ind1;
            nodetrisindex[indnode]++;
            ind0++;
        }
    }
    for(int ind1=0;ind1<(num_nodes);ind1++){
        nodetrisindex[ind1]=0;
    }
    for(int ind1=1;ind1<num_nodes;ind1++){
        nodetrisindex[ind1]=nodetrisindex[ind1-1]+num_trispernode[ind1-1];
    }

    //Calculate the number of connected edges
    num_rwgs=0;
    for(int ind1=0;ind1<num_tris;ind1++){
        for(int ind2=0;ind2<3;ind2++){
            if(ind2==0){
                nodes[0]=tri_nodes[ind1*3+0];
                nodes[1]=tri_nodes[ind1*3+1];
                nodes[2]=tri_nodes[ind1*3+2];

            } else if(ind2==1){
                nodes[0]=tri_nodes[ind1*3+1];
                nodes[1]=tri_nodes[ind1*3+2];
                nodes[2]=tri_nodes[ind1*3+0];
            } else if(ind2==2){
                nodes[0]=tri_nodes[ind1*3+2];
                nodes[1]=tri_nodes[ind1*3+0];
                nodes[2]=tri_nodes[ind1*3+1];
            }
            for(int indtri1=0;indtri1<num_trispernode[nodes[0]];indtri1++){
                tri1=nodetris[nodetrisindex[nodes[0]]+indtri1];
                for(int indtri2=0;indtri2<num_trispernode[nodes[1]];indtri2++){
                    tri2=nodetris[nodetrisindex[nodes[1]]+indtri2];
                    if((tri1==tri2) && (ind1!=tri1) && (ind1<tri1)){
                        num_rwgs++;
                    }
                }
            }
        }
    }
    printf("Number of Edges: %d \n", num_rwgs);
    rwg_tris = new int[num_rwgs*2];
    rwg_nodes = new int[num_rwgs*4];
    rwg_lengths = new double[num_rwgs];


    // Determine the common edges between triangles
    num_rwgs=0;
    for(int ind1=0;ind1<num_tris;ind1++){
        for(int ind2=0;ind2<3;ind2++){
            if(ind2==0){
                nodes[0]=tri_nodes[ind1*3+0];
                nodes[1]=tri_nodes[ind1*3+1];
                nodes[2]=tri_nodes[ind1*3+2];

            } else if(ind2==1){
                nodes[0]=tri_nodes[ind1*3+1];
                nodes[1]=tri_nodes[ind1*3+2];
                nodes[2]=tri_nodes[ind1*3+0];
            } else if(ind2==2){
                nodes[0]=tri_nodes[ind1*3+2];
                nodes[1]=tri_nodes[ind1*3+0];
                nodes[2]=tri_nodes[ind1*3+1];
            }
            for(int indtri1=0;indtri1<num_trispernode[nodes[0]];indtri1++){
                tri1=nodetris[nodetrisindex[nodes[0]]+indtri1];
                for(int indtri2=0;indtri2<num_trispernode[nodes[1]];indtri2++){
                    tri2=nodetris[nodetrisindex[nodes[1]]+indtri2];
                    if((tri1==tri2) && (ind1!=tri1) && (ind1<tri1)){

                        rwg_tris[num_rwgs*2+0]=ind1;
                        rwg_tris[num_rwgs*2+1]=tri1;
                        rwg_nodes[num_rwgs*4+0]=nodes[0];
                        rwg_nodes[num_rwgs*4+1]=nodes[1];
                        rwg_nodes[num_rwgs*4+2]=nodes[2];
                        element1node[0]=tri_nodes[tri1*3+0];
                        element1node[1]=tri_nodes[tri1*3+1];
                        element1node[2]=tri_nodes[tri1*3+2];
                        if((element1node[0]!=nodes[0]) && (element1node[0]!=nodes[1]) && (element1node[0]!=nodes[2]))
                            rwg_nodes[num_rwgs*4+3]=tri_nodes[tri1*3+0];
                        if((element1node[1]!=nodes[0]) && (element1node[1]!=nodes[1]) && (element1node[1]!=nodes[2]))
                            rwg_nodes[num_rwgs*4+3]=tri_nodes[tri1*3+1];
                        if((element1node[2]!=nodes[0]) && (element1node[2]!=nodes[1]) && (element1node[2]!=nodes[2]))
                            rwg_nodes[num_rwgs*4+3]=tri_nodes[tri1*3+2];

                        rwg_nodepos1[0]=node_pos[nodes[0]*3+0];
                        rwg_nodepos1[1]=node_pos[nodes[0]*3+1];
                        rwg_nodepos1[2]=node_pos[nodes[0]*3+2];
                        rwg_nodepos2[0]=node_pos[nodes[1]*3+0];
                        rwg_nodepos2[1]=node_pos[nodes[1]*3+1];
                        rwg_nodepos2[2]=node_pos[nodes[1]*3+2];
                        rwg_lengths[num_rwgs]=sqrt((rwg_nodepos1[0]-rwg_nodepos2[0])*(rwg_nodepos1[0]-rwg_nodepos2[0])\
                                                  +(rwg_nodepos1[1]-rwg_nodepos2[1])*(rwg_nodepos1[1]-rwg_nodepos2[1])\
                                                  +(rwg_nodepos1[2]-rwg_nodepos2[2])*(rwg_nodepos1[2]-rwg_nodepos2[2]));
                        num_rwgs++;
                    }
                }
            }
        }
    }
    delete[] nodetrisindex; nodetrisindex=nullptr;
    delete[] nodetris;
    delete[] num_trispernode;

    // Generate Triangle to RWG mapping
    tri2rwgindex = new int[num_tris*3];
    tri2rwgnode = new int[num_tris*3];
    tri2rwgsign = new int[num_tris*3];
    tri_areas = new double[num_tris];
    vec_trinormals = new double[num_tris*3];

    for(int ind1=0;ind1<(num_tris*3);ind1++){
        tri2rwgindex[ind1]=-1;
        tri2rwgnode[ind1]=-1;
        tri2rwgsign[ind1]=-1;
    }
    for(int ind1=0;ind1<num_rwgs;ind1++){
        tri1=rwg_tris[ind1*2+0];
        ind0=0;
        while(tri2rwgindex[tri1*3+ind0]!=-1)
            ind0++;
        tri2rwgindex[tri1*3+ind0]=ind1;
        tri2rwgnode[tri1*3+ind0]=rwg_nodes[ind1*4+2];
        tri2rwgsign[tri1*3+ind0]=1;
        tri1=rwg_tris[ind1*2+1];
        ind0=0;
        while(tri2rwgindex[tri1*3+ind0]!=-1)
            ind0++;
        tri2rwgindex[tri1*3+ind0]=ind1;
        tri2rwgnode[tri1*3+ind0]=rwg_nodes[ind1*4+3];
        tri2rwgsign[tri1*3+ind0]=-1;
    }

    // Calculate Triangle Areas and Surface Normals
    for(int ind1=0;ind1<num_tris;ind1++){
        nodes[0]=tri_nodes[ind1*3+0];
        nodes[1]=tri_nodes[ind1*3+1];
        nodes[2]=tri_nodes[ind1*3+2];
        temp_node_pos1[0]=node_pos[nodes[0]*3+0];
        temp_node_pos1[1]=node_pos[nodes[0]*3+1];
        temp_node_pos1[2]=node_pos[nodes[0]*3+2];
        temp_node_pos2[0]=node_pos[nodes[1]*3+0];
        temp_node_pos2[1]=node_pos[nodes[1]*3+1];
        temp_node_pos2[2]=node_pos[nodes[1]*3+2];
        temp_node_pos3[0]=node_pos[nodes[2]*3+0];
        temp_node_pos3[1]=node_pos[nodes[2]*3+1];
        temp_node_pos3[2]=node_pos[nodes[2]*3+2];
        temp_vec1[0]=temp_node_pos2[0]-temp_node_pos1[0];
        temp_vec1[1]=temp_node_pos2[1]-temp_node_pos1[1];
        temp_vec1[2]=temp_node_pos2[2]-temp_node_pos1[2];
        temp_vec2[0]=temp_node_pos3[0]-temp_node_pos1[0];
        temp_vec2[1]=temp_node_pos3[1]-temp_node_pos1[1];
        temp_vec2[2]=temp_node_pos3[2]-temp_node_pos1[2];
        // Cross-product
        temp_crossprod[0]=temp_vec1[1]*temp_vec2[2]-temp_vec1[2]*temp_vec2[1];
        temp_crossprod[1]=temp_vec1[2]*temp_vec2[0]-temp_vec1[0]*temp_vec2[2];
        temp_crossprod[2]=temp_vec1[0]*temp_vec2[1]-temp_vec1[1]*temp_vec2[0];
        temp_normalize=sqrt(temp_crossprod[0]*temp_crossprod[0]+temp_crossprod[1]*temp_crossprod[1]+temp_crossprod[2]*temp_crossprod[2]);
        tri_areas[ind1]=temp_normalize/2;
        vec_trinormals[ind1*3+0]=temp_crossprod[0]/temp_normalize;
        vec_trinormals[ind1*3+1]=temp_crossprod[1]/temp_normalize;
        vec_trinormals[ind1*3+2]=temp_crossprod[2]/temp_normalize;

    }

    // Check geometry for open edges
    for(int ind1=0;ind1<num_tris*3;ind1++){
        if(tri2rwgindex[ind1]==-1){
            opengeometry=1;
            printf("Open triangle is: %d \n", ind1/3);
            break;
        }
    }

    if(opengeometry==1)
        printf("Geometry is open! \n");
    else
        printf("Geometry is closed! \n");
    printf("\n");

    std::cout.flush();

}
