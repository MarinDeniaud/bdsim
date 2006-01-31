# ALL LENGTHS MUST BE GIVEN IN MILLIMETRES
# ALL ROTATIONS MUST BE GIVEN IN RADIANS

DROP DATABASE IF EXISTS INVOUTERMASKS;
CREATE DATABASE INVOUTERMASKS;
USE INVOUTERMASKS;

CREATE TABLE INVHCAL_CONE (
    POSX          DOUBLE(10,3), #
    POSY          DOUBLE(10,3), #
    POSZ          DOUBLE(10,3), #
    RED           DOUBLE(10,3), #
    GREEN         DOUBLE(10,3), #
    BLUE          DOUBLE(10,3), #
    VISATT        VARCHAR(32),  # I = INVISIBLE, S = SOLID, W = WIREFRAME
    LENGTH        DOUBLE(10,3), #
    RINNERSTART   DOUBLE(10,3), #
    RINNEREND     DOUBLE(10,3), #
    ROUTERSTART   DOUBLE(10,3), #
    ROUTEREND     DOUBLE(10,3), #
    ROTPSI        DOUBLE(10,3), #
    ROTTHETA      DOUBLE(10,3), #
    ROTPHI        DOUBLE(10,3), #
    MATERIAL      VARCHAR(32),  # MATERIAL, CGA LITERAL NAME
    NAME          VARCHAR(32)   # NAME OF SOLID, LOGICAL, AND PHYSICAL VOLUME
);
                      
INSERT INTO INVHCAL_CONE VALUES (3165e-3, 0.0, 8830.0, 0.0, 0.8, 0.0, "S", 330, 350.0, 350.0, 500.0, 500.0, 0., 0., -1e-3, "TUNGSTEN", ""); 

INSERT INTO INVHCAL_CONE VALUES (3790e-3, 0.0, 8205.0, 0.0, 0.8, 0.0, "S", 920, 300.0, 300.0, 500.0, 500.0, 0., 0., -1e-3, "TUNGSTEN", ""); 

#############################################################

CREATE TABLE INVPOLETIP_CONE (
    POSX          DOUBLE(10,3), #
    POSY          DOUBLE(10,3), #
    POSZ          DOUBLE(10,3), #
    RED           DOUBLE(10,3), #
    GREEN         DOUBLE(10,3), #
    BLUE          DOUBLE(10,3), #
    VISATT        VARCHAR(32),  # I = INVISIBLE, S = SOLID, W = WIREFRAME
    LENGTH        DOUBLE(10,3), #
    RINNERSTART   DOUBLE(10,3), #
    RINNEREND     DOUBLE(10,3), #
    ROUTERSTART   DOUBLE(10,3), #
    ROUTEREND     DOUBLE(10,3), #
    ROTPSI        DOUBLE(10,3), #
    ROTTHETA      DOUBLE(10,3), #
    ROTPHI        DOUBLE(10,3), #
    MATERIAL      VARCHAR(32),  # MATERIAL, CGA LITERAL NAME
    NAME          VARCHAR(32)   # NAME OF SOLID, LOGICAL, AND PHYSICAL VOLUME
);
                      
INSERT INTO INVPOLETIP_CONE VALUES (4626e-3, 0.0, 7370.0, 1.0, 0.0, 0.0, "S", 750, 450.0, 350.0, 500.0, 500.0, 0., 0., -1e-3, "IRON", ""); 

#############################################################

CREATE TABLE INVECAL_CONE (
    POSX          DOUBLE(10,3), #
    POSY          DOUBLE(10,3), #
    POSZ          DOUBLE(10,3), #
    RED           DOUBLE(10,3), #
    GREEN         DOUBLE(10,3), #
    BLUE          DOUBLE(10,3), #
    VISATT        VARCHAR(32),  # I = INVISIBLE, S = SOLID, W = WIREFRAME
    LENGTH        DOUBLE(10,3), #
    RINNERSTART   DOUBLE(10,3), #
    RINNEREND     DOUBLE(10,3), #
    ROUTERSTART   DOUBLE(10,3), #
    ROUTEREND     DOUBLE(10,3), #
    ROTPSI        DOUBLE(10,3), #
    ROTTHETA      DOUBLE(10,3), #
    ROTPHI        DOUBLE(10,3), #
    MATERIAL      VARCHAR(32),  # MATERIAL, CGA LITERAL NAME
    NAME          VARCHAR(32)   # NAME OF SOLID, LOGICAL, AND PHYSICAL VOLUME
);

INSERT INTO INVECAL_CONE VALUES (2900e-3, 0.0, 9095.0, 0.2, 0.6, 0.2, "S", 200, 250.0, 250.0, 500.0, 500.0, 0., 0., -1e-3, "TUNGSTEN", ""); 
