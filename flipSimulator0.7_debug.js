// flipSimulator.js propiedades grafico 50 step en 712 y 728 eliminado tangencial slip grafico en logaritmo

// --- Estructura de una Partícula (sin cambios, sigue igual) ---
/**
 * @typedef {Object} Particle
 * @property {number} x - Posición X en metros.
 * @property {number} y - Posición Y en metros.
 * @property {number} vx - Velocidad X en m/s.
 * @property {number} vy - Velocidad Y en m/s.
 */

// --- Clase para la Rejilla MAC (MACGrid) (sin cambios, sigue igual) ---
class MACGrid {
    constructor(nx, ny, cellSize) {
        this.nx = nx;
        this.ny = ny;
        this.cellSize = cellSize;

        this.u = Array(nx + 1).fill(null).map(() => Array(ny).fill(0.0));
        this.v = Array(nx).fill(null).map(() => Array(ny + 1).fill(0.0));
        this.p = Array(nx).fill(null).map(() => Array(ny).fill(0.0));
        this.isFluid = Array(nx).fill(null).map(() => Array(ny).fill(0));
        this.isSolid = Array(nx).fill(null).map(() => Array(ny).fill(0)); // 0=no sólido, 1=sólido

        this.particleCountPerCell = Array(nx).fill(null).map(() => Array(ny).fill(0));

        this.isEmitterCell = Array(nx).fill(null).map(() => Array(ny).fill(0));

        this.u_before_pressure = Array(nx + 1).fill(null).map(() => Array(ny).fill(0.0));
        this.v_before_pressure = Array(nx).fill(null).map(() => Array(ny + 1).fill(0.0));
        this.cellDensity = Array(nx).fill(null).map(() => Array(ny).fill(0.0));
        
    }
    // ... (otros métodos de MACGrid si los tuvieras)
}

// --- Clase FLIPSimulator (Aquí es donde insertamos el nuevo método) ---
class FLIPSimulator {
    constructor(domainWidthMeters, domainHeightMeters, cellSizeMeters) {
        this.GRAVITY_Y = -9.81; 
        this.SIM_DENSITY = 1.0; 

        this.particles = []; 
        const nx = Math.floor(domainWidthMeters / cellSizeMeters);
        const ny = Math.floor(domainHeightMeters / cellSizeMeters);
        this.grid = new MACGrid(nx, ny, cellSizeMeters);

        this.picFlipAlpha = 0.95; 
        this.restitutionCoefficient = 0.3; 
        this.frictionFactor = 0.95; 
        this.numPressureIterations = 8000; // Lo ajustas en el HTML
        this.pressureRelaxationFactor = 1; 
        this.pressureTolerance = 1e-6;
        this.useNoSlip = false; // Para decidir si aplicar no-slip tangencial
        
        // --- NUEVAS/MODIFICADAS PROPIEDADES PARA SIEMBRA DE PARTÍCULAS ---
        this.particlesPerCellDim = 2; // Para NxN partículas, ej. 2 significa 2x2=4 partículas
        this.particleJitterStrength = 0.25; // Fuerza del jitter (0 a 1, donde 0.5 es la mitad del espaciado)
        // --- FIN NUEVAS/MODIFICADAS ---

            // NUEVO: Parámetros para el control de densidad
        // Si siembras 2x2=4 partículas por celda, un buen punto de partida para la densidad objetivo
        // podría ser cercano a la suma de pesos que estas 4 partículas aportarían si están bien centradas.
        // El peso de una partícula perfectamente centrada es 1. Si las 4 están distribuidas,
        // la suma de pesos puede ser menor que 4. Experimenta con este valor.
        // Comencemos con un valor que represente aproximadamente 2-3 partículas "efectivas".
        //this.targetCellDensity = 5.0; // Con celdas de 2x2 ha funcionado bien 5
        const numParticlesNominalPerCell = this.particlesPerCellDim * this.particlesPerCellDim;
        const targetDensityFactor = (0.25 * this.particlesPerCellDim) + 0.75;

        this.targetCellDensity = numParticlesNominalPerCell * targetDensityFactor; // o numParticlesNominalPerCell * 1.25
        this.densityCorrectionFactor = 0.1; // VALOR PARA EXPERIMENTAR (como en el 3D)
        //this.MAX_REASONABLE_CELL_DENSITY = 10.0;  // O 15.0, el valor del clamp que establecimos
        // MAX_REASONABLE_CELL_DENSITY debe ser mayor que targetCellDensity.
        // Por ejemplo, 2 veces el objetivo, o el objetivo + un margen, con un mínimo.
        /////// this.MAX_REASONABLE_CELL_DENSITY = Math.max(this.targetCellDensity + 5.0, this.targetCellDensity * 1.5);
        this.MAX_REASONABLE_CELL_DENSITY = 10.0;

        // Para targetCellDensity = 5.0 (con N=2):
        // Math.max(5.0 + 5.0, 5.0 * 1.5) = Math.max(10.0, 7.5) = 10.0.
        // Si N=1, numParticlesNominal=1, targetCellDensity=1.25.
        // MAX_REASONABLE_CELL_DENSITY = Math.max(1.25 + 5.0, 1.25 * 1.5) = Math.max(6.25, 1.875) = 6.25.
        // Si N=3, numParticlesNominal=9, targetCellDensity=11.25.
        // MAX_REASONABLE_CELL_DENSITY = Math.max(11.25 + 5.0, 11.25 * 1.5) = Math.max(16.25, 16.875) = 16.875.
        // Esto parece un escalado razonable. Asegúrate que el clamp en P2G usa este valor.

        //console.log("DEBUG CONSTRUCTOR - Initial this.particlesPerCellDim:", this.particlesPerCellDim);
        //console.log("DEBUG CONSTRUCTOR - numParticlesNominalPerCell:", numParticlesNominalPerCell);
        //console.log("DEBUG CONSTRUCTOR - Calculated this.targetCellDensity:", this.targetCellDensity);
        //console.log("DEBUG CONSTRUCTOR - Calculated this.MAX_REASONABLE_CELL_DENSITY:", this.MAX_REASONABLE_CELL_DENSITY);
        //console.log("DEBUG CONSTRUCTOR - this.densityCorrectionFactor:", this.densityCorrectionFactor);

        this.solidAABBs = []; // Array para almacenar los AABBs de los sólidos {xMin, xMax, yMin, yMax}
        this.hasComplexSolids = false; // Por defecto, no hay sólidos que requieran sub-stepping

        this.pressureIterationsData = [];

        this.simulationTime = 0.0;
        this.maxHeightData = []; 
        this.particleCountHistoryData = []; // NUEVO ARRAY
        this.graphTimeWindowSeconds = 80;

        // --- NUEVO: Banderas para controlar contornos abiertos ---
        this.boundaryOpen = {
            left:   false, // Poner a true para que la pared izquierda sea abierta
            right:  false, // Poner a true para que la pared derecha sea abierta
            top:    false, // Poner a true para que la pared superior sea abierta
            bottom: false  // Poner a true para que la pared inferior sea abierta (cuidado con la gravedad)
        };
        // --- FIN NUEVO ---


        console.log(`Simulador inicializado: ${nx}x${ny} celdas, tamaño de celda: ${cellSizeMeters.toFixed(3)}m`);
        console.log(`Dominio físico: ${this.getDomainWidthMeters()}m x ${this.getDomainHeightMeters()}m`);
    }

    getDomainWidthMeters() {
        return this.grid.nx * this.grid.cellSize;
    }

    getDomainHeightMeters() {
        return this.grid.ny * this.grid.cellSize;
    }

    addParticle(x_m, y_m, vx_mps = 0, vy_mps = 0) {
        this.particles.push({ x: x_m, y: y_m, vx: vx_mps, vy: vy_mps });
    }
    
    addRotatedSolidBox(centerX_m, centerY_m, width_m, height_m, angle_degrees) {
        const h = this.grid.cellSize;
        const nx = this.grid.nx;
        const ny = this.grid.ny;

        this.hasComplexSolids = true; // Indicar que ahora hay geometría que podría necesitar sub-pasos


        const angle_rad = angle_degrees * (Math.PI / 180.0);
        const cos_a = Math.cos(angle_rad);
        const sin_a = Math.sin(angle_rad);

        // Para optimizar, solo iteraremos sobre las celdas que podrían
        // estar cubiertas por el rectángulo rotado. Calculamos un AABB
        // que envuelva al rectángulo rotado.
        const half_w = width_m / 2.0;
        const half_h = height_m / 2.0;

        // Vértices del rectángulo NO rotado, relativos a su centro (0,0)
        const corners = [
            { x: -half_w, y: -half_h },
            { x:  half_w, y: -half_h },
            { x:  half_w, y:  half_h },
            { x: -half_w, y:  half_h }
        ];

        let minX_rot = Infinity, maxX_rot = -Infinity;
        let minY_rot = Infinity, maxY_rot = -Infinity;

        // Rotar los vértices y encontrar el AABB del rectángulo rotado
        for (const corner of corners) {
            const rotatedX = corner.x * cos_a - corner.y * sin_a + centerX_m;
            const rotatedY = corner.x * sin_a + corner.y * cos_a + centerY_m;
            minX_rot = Math.min(minX_rot, rotatedX);
            maxX_rot = Math.max(maxX_rot, rotatedX);
            minY_rot = Math.min(minY_rot, rotatedY);
            maxY_rot = Math.max(maxY_rot, rotatedY);
        }

        // Convertir el AABB del rotado a índices de celda
        const iStart = Math.max(0, Math.floor(minX_rot / h));
        const iEnd = Math.min(nx - 1, Math.ceil(maxX_rot / h)); // Usar ceil para asegurar que cubrimos todo
        const jStart = Math.max(0, Math.floor(minY_rot / h));
        const jEnd = Math.min(ny - 1, Math.ceil(maxY_rot / h));

        // Iterar sobre las celdas dentro de este AABB
        for (let i = iStart; i <= iEnd; i++) {
            for (let j = jStart; j <= jEnd; j++) {
                // Coordenadas del centro de la celda actual P[i][j]
                const cellCenterX_m = (i + 0.5) * h;
                const cellCenterY_m = (j + 0.5) * h;

                // Para verificar si el centro de la celda está dentro del rectángulo ROTADO,
                // es más fácil rotar el centro de la celda EN SENTIDO CONTRARIO
                // y ver si cae dentro del rectángulo original NO ROTADO (que está centrado en el origen).
                
                // 1. Trasladar el centro de la celda para que el centro del rectángulo sea el origen
                const translatedX = cellCenterX_m - centerX_m;
                const translatedY = cellCenterY_m - centerY_m;

                // 2. Rotar el centro de la celda trasladado por -angle_rad
                const localX = translatedX * cos_a + translatedY * sin_a; // cos(-a)=cos(a), sin(-a)=-sin(a)
                const localY = -translatedX * sin_a + translatedY * cos_a;

                // 3. Comprobar si el punto (localX, localY) está dentro del rectángulo NO rotado
                //    (que está centrado en el origen y tiene dimensiones width_m, height_m)
                if (Math.abs(localX) <= half_w + 1e-6 && Math.abs(localY) <= half_h + 1e-6) { // Añadido un pequeño epsilon por precisión
                    if (i >= 0 && i < nx && j >= 0 && j < ny) { // Asegurar que el índice de celda es válido
                        this.grid.isSolid[i][j] = 1;
                        this.grid.isFluid[i][j] = 0; 
                    }
                }
            }
        }
        console.log(`Caja sólida ROTADA añadida (rasterizada en celdas).`);
        // Nota: this.solidAABBs NO se actualiza con un AABB para esta forma rotada,
        // ya que la colisión de partículas actual en advectParticles se basa en grid.isSolid.
        // Si quisiéramos una colisión más precisa con el AABB del rotado, necesitaríamos
        // pasar esa información y modificar advectParticles.
    }
    
    // ... (resto de los métodos: particleToGrid, applyExternalForces, etc.) ...
    // Asegúrate de que tu método advectParticles sea la versión que comprueba colisiones
    // con this.grid.isSolid[ci][cj] para los obstáculos internos.


    addSolidBox(xMin_m, yMin_m, width_m, height_m) {
        const h = this.grid.cellSize;
        const nx = this.grid.nx;
        const ny = this.grid.ny;
    
        const iStart = Math.max(0, Math.floor(xMin_m / h));
        const iEnd = Math.min(nx - 1, Math.floor((xMin_m + width_m - 1e-9) / h)); // Usar 1e-9 para ser más preciso con el borde
        const jStart = Math.max(0, Math.floor(yMin_m / h));
        const jEnd = Math.min(ny - 1, Math.floor((yMin_m + height_m - 1e-9) / h));
    
        for (let i = iStart; i <= iEnd; i++) {
            for (let j = jStart; j <= jEnd; j++) {
                if (i >= 0 && i < nx && j >= 0 && j < ny) {
                    this.grid.isSolid[i][j] = 1;
                    this.grid.isFluid[i][j] = 0; 
                }
            }
        }
        
        // Almacenar las coordenadas exactas del AABB del obstáculo en metros
        const boxMinX = iStart * h;
        const boxMaxX = (iEnd + 1) * h;
        const boxMinY = jStart * h;
        const boxMaxY = (jEnd + 1) * h;
        this.solidAABBs.push({ 
            xMin: boxMinX, xMax: boxMaxX, 
            yMin: boxMinY, yMax: boxMaxY 
        });
        console.log(`Caja sólida añadida de celda (${iStart},${jStart}) a (${iEnd},${jEnd})`);
        console.log(`AABB Sólido: x[${boxMinX.toFixed(3)}-${boxMaxX.toFixed(3)}], y[${boxMinY.toFixed(3)}-${boxMaxY.toFixed(3)}]`);
    }

    addSolidDiagonal(x1_m, y1_m, x2_m, y2_m, thickness_cells = 10) {
        const h = this.grid.cellSize;
        const nx = this.grid.nx;
        const ny = this.grid.ny;
    
        // Convertir a coordenadas de celda
        let c1x = Math.floor(x1_m / h);
        let c1y = Math.floor(y1_m / h);
        let c2x = Math.floor(x2_m / h);
        let c2y = Math.floor(y2_m / h);
    
        // Algoritmo de línea de Bresenham (simplificado para marcar celdas)
        let dx = Math.abs(c2x - c1x);
        let dy = -Math.abs(c2y - c1y);
        let sx = c1x < c2x ? 1 : -1;
        let sy = c1y < c2y ? 1 : -1;
        let err = dx + dy; 
        let e2;
    
        console.log(`Dibujando línea sólida de (${c1x},${c1y}) a (${c2x},${c2y})`);
    
        while (true) {
            // Marcar la celda actual y su vecindad según thickness_cells
            for (let thx = -Math.floor((thickness_cells-1)/2); thx <= Math.ceil((thickness_cells-1)/2); thx++) {
                for (let thy = -Math.floor((thickness_cells-1)/2); thy <= Math.ceil((thickness_cells-1)/2); thy++) {
                    const curX = c1x + thx;
                    const curY = c1y + thy;
                    if (curX >= 0 && curX < nx && curY >= 0 && curY < ny) {
                        this.grid.isSolid[curX][curY] = 1;
                        this.grid.isFluid[curX][curY] = 0;
                    }
                }
            }
            
            if (c1x === c2x && c1y === c2y) break;
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; c1x += sx; }
            if (e2 <= dx) { err += dx; c1y += sy; }
        }
        console.log("Línea sólida diagonal aproximada añadida.");
    }

    initializeFluidVolume(xMin_m, yMin_m, width_m, height_m) { // Se eliminó particlesPerCellDim de los parámetros
        const xMax_m = xMin_m + width_m;
        const yMax_m = yMin_m + height_m;
        // const jitterStrength = 0.25; 

        for (let i = 0; i < this.grid.nx; i++) {
        for (let j = 0; j < this.grid.ny; j++) {
            const cellMinX = i * this.grid.cellSize;
            const cellMinY = j * this.grid.cellSize;
            const cellCenterX = cellMinX + this.grid.cellSize / 2;
            const cellCenterY = cellMinY + this.grid.cellSize / 2;

            if (cellCenterX >= xMin_m && cellCenterX < xMax_m &&
                cellCenterY >= yMin_m && cellCenterY < yMax_m) {

                // Comprobar si la celda actual (i,j) es sólida ANTES de intentar sembrar partículas
                if (this.grid.isSolid[i][j]) {
                    continue; // Saltar esta celda si es sólida
                }

                for (let px_idx = 0; px_idx < this.particlesPerCellDim; px_idx++) {
                    for (let py_idx = 0; py_idx < this.particlesPerCellDim; py_idx++) {
                        const spacing = this.grid.cellSize / this.particlesPerCellDim;
                        const jitterX = (Math.random() - 0.5) * spacing * this.particleJitterStrength;
                        const jitterY = (Math.random() - 0.5) * spacing * this.particleJitterStrength;
                        const particleX = cellMinX + (px_idx + 0.5) * spacing + jitterX;
                        const particleY = cellMinY + (py_idx + 0.5) * spacing + jitterY;

                        if (particleX >= xMin_m && particleX < xMax_m && 
                            particleY >= yMin_m && particleY < yMax_m) {
                            // No es necesario volver a comprobar si la celda es sólida aquí si ya lo hicimos arriba
                            this.addParticle(particleX, particleY, 0, 0);
                        }
                    }
                }
            }
        }
    }
        // console.log(`Volumen de fluido inicializado con ${this.particles.length} partículas.`);
    }

    addParticlesInCell(ci, cj, isSimRunning, emitVx = 0, emitVy = 0) { 
    if (ci < 0 || ci >= this.grid.nx || cj < 0 || cj >= this.grid.ny) {
        return false; 
    }
    if (this.grid.isSolid[ci][cj]) {
        return false; 
    }

    const h = this.grid.cellSize;
    if (!isSimRunning) {
        const cellMinX_clear = ci * h;
        const cellMaxX_clear = (ci + 1) * h;
        const cellMinY_clear = cj * h;
        const cellMaxY_clear = (cj + 1) * h;
        this.particles = this.particles.filter(p => {
            return !(p.x >= cellMinX_clear && p.x < cellMaxX_clear &&
                       p.y >= cellMinY_clear && p.y < cellMaxY_clear);
        });
    }

    let particlesActuallyAddedThisCall = 0;
    const cellMinX_seed = ci * h;
    const cellMinY_seed = cj * h;
    const expectedToAdd = this.particlesPerCellDim * this.particlesPerCellDim;

    for (let px_idx = 0; px_idx < this.particlesPerCellDim; px_idx++) {
        for (let py_idx = 0; py_idx < this.particlesPerCellDim; py_idx++) {
            const spacing = h / this.particlesPerCellDim;
            const jitterX = (Math.random() - 0.5) * spacing * this.particleJitterStrength;
            const jitterY = (Math.random() - 0.5) * spacing * this.particleJitterStrength;
            
            let particleX = cellMinX_seed + (px_idx + 0.5) * spacing + jitterX;
            let particleY = cellMinY_seed + (py_idx + 0.5) * spacing + jitterY;

            const safetyMargin = h * 0.01;
            particleX = Math.max(safetyMargin, Math.min(particleX, this.getDomainWidthMeters() - safetyMargin));
            particleY = Math.max(safetyMargin, Math.min(particleY, this.getDomainHeightMeters() - safetyMargin));
            
            const final_ci = Math.floor(particleX / h);
            const final_cj = Math.floor(particleY / h);

            if (final_ci === ci && final_cj === cj) { 
               // ***** MODIFICADO: Usar emitVx y emitVy *****
               this.addParticle(particleX, particleY, emitVx, emitVy);
               particlesActuallyAddedThisCall++;
            }
        }
    }

    // ... (tus logs de depuración si los mantienes) ...
    if (isSimRunning) { 
        if (particlesActuallyAddedThisCall > 0) {
            if (particlesActuallyAddedThisCall < expectedToAdd) {
                console.warn(`LOG DE ADD PARTICLES: Emisor/Adición en (${ci},${cj}) añadió ${particlesActuallyAddedThisCall} de ${expectedToAdd} esperadas. (t=${this.simulationTime.toFixed(3)})`);
            } else {
                // console.log(`LOG DE ADD PARTICLES: Emisor/Adición en (${ci},${cj}) añadió ${particlesActuallyAddedThisCall} partículas. (t=${this.simulationTime.toFixed(3)})`);
            }
        } else if (expectedToAdd > 0) { 
            // console.log(`LOG DE ADD PARTICLES: Emisor/Adición en (${ci},${cj}) NO añadió partículas (0 de ${expectedToAdd}). (t=${this.simulationTime.toFixed(3)})`);
        }
    }
    return particlesActuallyAddedThisCall > 0; 
}

    // Dentro de la clase FLIPSimulator:
    particleToGrid() {
        const particles = this.particles;
        const grid = this.grid;
        const nx = grid.nx;
        const ny = grid.ny;
        const h = grid.cellSize;
        const isSolid = grid.isSolid;
        //const MAX_REASONABLE_CELL_DENSITY = 10.0; // VALOR PARA EXPERIMENTAR. 
        
        console.log(`Número total de partículas (this.particles.length): ${this.particles.length}`);


        // 1. Crear e inicializar acumuladores temporales para este paso P2G
        const u_sum_vx = Array(nx + 1).fill(null).map(() => Array(ny).fill(0.0));
        const u_sum_weights = Array(nx + 1).fill(null).map(() => Array(ny).fill(0.0));
        const v_sum_vy = Array(nx).fill(null).map(() => Array(ny + 1).fill(0.0));
        const v_sum_weights = Array(nx).fill(null).map(() => Array(ny + 1).fill(0.0));
        const cell_sum_mass = Array(nx).fill(null).map(() => Array(ny).fill(0.0));

        // 2. Resetear velocidades de la rejilla y isFluid map
        for (let i = 0; i <= nx; i++) { for (let j = 0; j < ny; j++) { grid.u[i][j] = 0.0; } }
        for (let i = 0; i < nx; i++) { for (let j = 0; j <= ny; j++) { grid.v[i][j] = 0.0; } }
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                if (!isSolid[i][j]) {
                    grid.isFluid[i][j] = 0;
                } else {
                    grid.isFluid[i][j] = 0; // Sólida no es fluida y no debería marcarse aquí tampoco
                }
                grid.cellDensity[i][j] = 0.0; // NUEVO: Resetear cellDensity
                grid.particleCountPerCell[i][j] = 0; // NUEVO: Resetear contador de partículas

            }
        }

        // 3. Iterar sobre cada partícula para "esparcir" (scatter) su influencia
        for (const particle of particles) {
            const px = particle.x; const py = particle.y; const pvx = particle.vx; const pvy = particle.vy; const particleMass = 1.0;
            
            // --- A. Transferir particle.vx a las caras U cercanas ---
            const gux = px/h; const guy = py/h - 0.5; const i0u = Math.floor(gux); const j0u = Math.floor(guy); const fxu = gux - i0u; const fyu = guy - j0u;
            for (let jo = 0; jo <= 1; jo++) { 
                for (let io = 0; io <= 1; io++) { 
                    const ui = i0u + io; const uj = j0u + jo; 
                    if (ui >= 0 && ui <= nx && uj >= 0 && uj < ny) { 
                        const wux = (io === 0) ? (1 - fxu) : fxu; 
                        const wuy = (jo === 0) ? (1 - fyu) : fyu; 
                        const wu = wux * wuy; 
                        u_sum_vx[ui][uj] += pvx * wu; 
                        u_sum_weights[ui][uj] += wu; 
                    }
                }
            }

            // --- B. Transferir particle.vy a las caras V cercanas ---
            const gvx = px/h - 0.5; const gvy = py/h; const i0v = Math.floor(gvx); const j0v = Math.floor(gvy); const fxv = gvx - i0v; const fyv = gvy - j0v;
            for (let jo = 0; jo <= 1; jo++) { 
                for (let io = 0; io <= 1; io++) { 
                    const vi = i0v + io; const vj = j0v + jo; 
                    if (vi >= 0 && vi < nx && vj >= 0 && vj <= ny) { 
                        const wvx = (io === 0) ? (1 - fxv) : fxv; 
                        const wvy = (jo === 0) ? (1 - fyv) : fyv; 
                        const wv = wvx * wvy; 
                        v_sum_vy[vi][vj] += pvy * wv; 
                        v_sum_weights[vi][vj] += wv; 
                    }
                }
            }
            
            // --- C. Transferir "masa" a centros de celda P (para isFluid) ---
            const gpcx = px/h - 0.5; const gpcy = py/h - 0.5; 
            const i0c = Math.floor(gpcx); 
            const j0c = Math.floor(gpcy); 
            const fxc = gpcx - i0c;       
            const fyc = gpcy - j0c;       
            for (let j_offset = 0; j_offset <= 1; j_offset++) { 
                for (let i_offset = 0; i_offset <= 1; i_offset++) { 
                    const ci = i0c + i_offset; 
                    const cj = j0c + j_offset; 
                    if (ci >= 0 && ci < nx && cj >= 0 && cj < ny) { 
                        if (!isSolid[ci][cj]) { 
                            const wc_x = (i_offset === 0) ? (1 - fxc) : fxc; 
                            const wc_y = (j_offset === 0) ? (1 - fyc) : fyc; 
                            const weight_c = wc_x * wc_y; 
                            //cell_sum_mass[ci][cj] += particleMass * weight_c; 
                            grid.cellDensity[ci][cj] += particleMass * weight_c; // NUEVO: Acumular en grid.cellDensity
                            grid.cellDensity[ci][cj] = Math.min(grid.cellDensity[ci][cj], this.MAX_REASONABLE_CELL_DENSITY);

                        } 
                    } 
                }
            }
            // << --- LOS BUCLES DE NORMALIZACIÓN Y MARCADO DE isFluid IBAN AQUÍ INCORRECTAMENTE --- >>
         // --- FIN DEL BUCLE "for (const particle of particles)" ---

        // NUEVO: Contar en qué celda principal está la partícula
        const mainCi = Math.floor(px / h);
        const mainCj = Math.floor(py / h);
        if (mainCi >= 0 && mainCi < nx && mainCj >= 0 && mainCj < ny) {
            grid.particleCountPerCell[mainCi][mainCj]++;
        }
        }

        // --- 4. Normalizar velocidades en la rejilla (DESPUÉS de procesar TODAS las partículas) ---
        for (let i = 0; i <= nx; i++) { 
            for (let j = 0; j < ny; j++) { 
                if (u_sum_weights[i][j] > 1e-9) {
                    grid.u[i][j] = u_sum_vx[i][j] / u_sum_weights[i][j];
                } else {
                    grid.u[i][j] = 0.0; 
                }
            }
        }
        for (let i = 0; i < nx; i++) { 
            for (let j = 0; j <= ny; j++) { 
                if (v_sum_weights[i][j] > 1e-9) {
                    grid.v[i][j] = v_sum_vy[i][j] / v_sum_weights[i][j];
                } else {
                    grid.v[i][j] = 0.0; 
                }
            }
        }
        
        // --- 5. Marcar celdas fluidas (DESPUÉS de procesar TODAS las partículas) ---
        // Este bucle es el que realmente establece el estado final de isFluid para este timestep.
        for (let i = 0; i < nx; i++) { 
            for (let j = 0; j < ny; j++) { 
                if (!isSolid[i][j]) { // Solo considerar si la celda no es un sólido predefinido
                    //if (cell_sum_mass[i][j] > 0.01) { // Umbral para considerar una celda como fluida
                    if (grid.cellDensity[i][j] > 0.01) { // Umbral para considerar una celda como fluida
                        grid.isFluid[i][j] = 1; 
                    } else { 
                        grid.isFluid[i][j] = 0; 
                    } 
                } else { 
                    grid.isFluid[i][j] = 0; // Si es un sólido predefinido, no es fluido
                } 
            } 
        }
        /*
        // console.log("Módulo: P2G Completado.");
        if (this.simulationTime > 10 && this.simulationTime < 10.1) { // Imprimir solo en un intervalo corto para no inundar la consola
        let floorDensities = [];
        let floorParticleCounts = []; // Para comparar

        for (let i = 0; i < grid.nx; i++) {
            if (grid.isFluid[i][0]) { // Solo celdas fluidas del suelo
                floorDensities.push(grid.cellDensity[i][0].toFixed(3));
                floorParticleCounts.push(grid.particleCountPerCell[i][0]); // Añadir conteo
            }
        }
        if(floorDensities.length > 0) {
            console.log(`T=${this.simulationTime.toFixed(2)}s. SUELO Densidades: ${floorDensities.join(', ')}`);
            console.log(`T=${this.simulationTime.toFixed(2)}s. SUELO Counts:   ${floorParticleCounts.join(', ')}`);
        }
    }*/

    // --- NUEVO: Inspección de partículas para celdas con densidad anómala ---
    // --- NUEVO: Inspección de partículas para celdas con densidad anómala ---
    const ANOMALY_DENSITY_THRESHOLD = 50.0; 
    let anomaliesFoundThisFrame = 0; // Renombrado para evitar conflicto si 'anomaliesFound' es global
    const MAX_ANOMALIES_TO_LOG_PER_FRAME = 15;

    // Asegúrate de que nx, ny, h, grid, particles estén accesibles aquí.
    // Si este código está al final de particleToGrid, ya lo están.

    if (this.simulationTime > 10 && this.simulationTime < 10.1) { 
        for (let i_cell = 0; i_cell < nx; i_cell++) { // Usar nombres diferentes para los iteradores de celda
            for (let j_cell = 0; j_cell < ny; j_cell++) { 
                if (grid.cellDensity[i_cell][j_cell] > ANOMALY_DENSITY_THRESHOLD && anomaliesFoundThisFrame < MAX_ANOMALIES_TO_LOG_PER_FRAME) {
                    anomaliesFoundThisFrame++;
                    console.warn(`ANOMALÍA DE DENSIDAD en celda (${i_cell},${j_cell}): ${grid.cellDensity[i_cell][j_cell].toFixed(2)}. Partículas cercanas:`);
                    
                    let contributingParticlesCount = 0;
                    const MAX_PARTICLES_TO_LOG_PER_CELL = 20;

                    for (const currentParticle of particles) { // Usar 'currentParticle' para evitar confusión con 'particle' del bucle principal de P2G si existiera
                        const pX = currentParticle.x; // Definir pX y pY para la partícula actual de este bucle
                        const pY = currentParticle.y;

                        const gpcx = pX/h - 0.5; 
                        const gpcy = pY/h - 0.5;
                        const i0c = Math.floor(gpcx); 
                        const j0c = Math.floor(gpcy);

                        // Comprobar si la celda anómala (i_cell, j_cell) es una de las 4 celdas a las que la partícula contribuye
                        if ((i0c === i_cell && j0c === j_cell) || (i0c + 1 === i_cell && j0c === j_cell) ||
                            (i0c === i_cell && j0c + 1 === j_cell) || (i0c + 1 === i_cell && j0c + 1 === j_cell)) {
                            
                            if (contributingParticlesCount < MAX_PARTICLES_TO_LOG_PER_CELL) {
                                console.log(`  P(${contributingParticlesCount}): x=${pX.toFixed(3)}, y=${pY.toFixed(3)}`);
                            }
                            contributingParticlesCount++;
                        }
                    }
                    if (contributingParticlesCount >= MAX_PARTICLES_TO_LOG_PER_CELL) {
                        console.log(`  ... y ${contributingParticlesCount - MAX_PARTICLES_TO_LOG_PER_CELL} más partículas contribuyendo a celda (${i_cell},${j_cell}).`);
                    }
                    if (contributingParticlesCount === 0 && grid.cellDensity[i_cell][j_cell] > 0) { // Si la densidad es >0 pero no encontramos partículas, la lógica de detección de contribuyentes podría tener un caso borde
                         console.log(`  Advertencia: Densidad > 0 pero no se listaron partículas contribuyentes para celda (${i_cell},${j_cell}). Revisar lógica de detección de stencil.`);
                    }
                }
                if (anomaliesFoundThisFrame >= MAX_ANOMALIES_TO_LOG_PER_FRAME) break; 
            }
            if (anomaliesFoundThisFrame >= MAX_ANOMALIES_TO_LOG_PER_FRAME) break; 
        }
    
    }

    }

    enforceSolidCellBoundaryConditions() {
        const grid = this.grid;
        const u = grid.u;
        const v = grid.v;
        const isSolid = grid.isSolid;
        const nx = grid.nx;
        const ny = grid.ny;

        // Iterar sobre todas las celdas de la rejilla principal
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                if (isSolid[i][j]) {
                    // Esta celda (i,j) es sólida.
                    // Forzar a cero las velocidades en sus 4 caras.

                    // Cara-u izquierda de la celda sólida (i,j)
                    // Esta cara también es la cara-u derecha de la celda (i-1,j) si i > 0
                    if (i < nx) { // u[i][j] es válida hasta i = nx-1 para la celda (i-1,j)
                                // y hasta i=nx para la celda (i,j) como cara izquierda
                        u[i][j] = 0.0;
                    }
                    

                    // Cara-u derecha de la celda sólida (i,j)
                    // Esta cara es u[i+1][j]
                    if (i + 1 <= nx) {
                        u[i+1][j] = 0.0;
                    }

                    // Cara-v inferior de la celda sólida (i,j)
                    // Esta cara es v[i][j]
                    if (j < ny) { // v[i][j] es válida hasta j = ny-1 para la celda (i,j-1)
                                // y hasta j=ny para la celda (i,j) como cara inferior
                        v[i][j] = 0.0;
                    }

                    // Cara-v superior de la celda sólida (i,j)
                    // Esta cara es v[i][j+1]
                    if (j + 1 <= ny) {
                        v[i][j+1] = 0.0;
                    }
                }
            }
        }
        // console.log("Módulo: Condiciones de Contorno en Celdas Sólidas Aplicadas.");
    }

    applyExternalForces(dt) {
        // ... (Código completo de applyExternalForces como lo teníamos) ...
        const grid = this.grid; const gravityY = this.GRAVITY_Y;
        for (let i = 0; i < grid.nx; i++) { for (let j = 0; j <= grid.ny; j++) { grid.v[i][j] += gravityY * dt; }}
    }
/*
    enforceMACBoundaryConditions() { // Solo normales
        // ... (Código completo de enforceMACBoundaryConditions (solo normales) como lo teníamos) ...
        const grid = this.grid; const u = grid.u; const v = grid.v; const nx = grid.nx; const ny = grid.ny;
        for (let j = 0; j < ny; j++) { u[0][j] = 0.0; u[nx][j] = 0.0; }
        for (let i = 0; i < nx; i++) { v[i][0] = 0.0; v[i][ny] = 0.0; }
    }
*/

    enforceMACBoundaryConditions() {
        const grid = this.grid;
        const u = grid.u;
        const v = grid.v;
        const nx = grid.nx;
        const ny = grid.ny;

        // --- PAREDES VERTICALES (afectan componente normal u) ---
        // u[i][j] es la cara entre P[i-1][j] y P[i][j]
        // u[0][j] es la cara más a la izquierda. u[nx][j] es la cara más a la derecha.

        for (let j = 0; j < ny; j++) { // Itera sobre todas las alturas de las caras u
            // Límite Izquierdo (u[0][j])
            if (this.boundaryOpen.left) {
                if (nx > 0) u[0][j] = u[1][j]; // Condición de gradiente cero (outflow)
                // else u[0][j] = 0; // Caso improbable para una rejilla válida
            } else {
                u[0][j] = 0.0; // Pared cerrada
            }

            // Límite Derecho (u[nx][j])
            if (this.boundaryOpen.right) {
                if (nx > 0) u[nx][j] = u[nx - 1][j]; // Condición de gradiente cero (outflow)
                // else u[nx][j] = 0; // Caso improbable
            } else {
                u[nx][j] = 0.0; // Pared cerrada
            }
        }

        // --- PAREDES HORIZONTALES (afectan componente normal v) ---
        // v[i][j] es la cara entre P[i][j-1] y P[i][j]
        // v[i][0] es la cara más inferior. v[i][ny] es la cara más superior.

        for (let i = 0; i < nx; i++) { // Itera sobre todas las anchuras de las caras v
            // Límite Inferior (v[i][0])
            if (this.boundaryOpen.bottom) {
                if (ny > 0) v[i][0] = v[i][1]; // Condición de gradiente cero (outflow)
                // else v[i][0] = 0; // Caso improbable
            } else {
                v[i][0] = 0.0; // Pared cerrada
            }

            // Límite Superior (v[i][ny])
            if (this.boundaryOpen.top) {
                if (ny > 0) v[i][ny] = v[i][ny - 1]; // Condición de gradiente cero (outflow)
                // else v[i][ny] = 0; // Caso improbable
            } else {
                v[i][ny] = 0.0; // Pared cerrada
            }
        }
    }
 /*
    applyTangentialNoSlip() {
        const grid = this.grid;
        const u = grid.u;
        const v = grid.v;
        const isSolid = grid.isSolid; // Necesitamos isSolid
        const nx = grid.nx;
        const ny = grid.ny;

        if (!this.useNoSlip) { // Si useNoSlip es false, no hacer nada
            return;
        }

        // --- NO-SLIP TANGENCIAL EN LÍMITES DEL DOMINIO ---
        // Paredes horizontales del dominio (superior e inferior) - afecta a 'u'
        for (let i = 0; i <= nx; i++) {
            u[i][0] = 0.0;    // Tangencial en la pared inferior j=0
            if (ny > 0) {
                u[i][ny-1] = 0.0; // Tangencial en la pared superior j=ny (afecta a la cara u más cercana)
            }
        }
        // Paredes verticales del dominio (izquierda y derecha) - afecta a 'v'
        for (let j = 0; j <= ny; j++) {
            v[0][j] = 0.0;    // Tangencial en la pared izquierda i=0
            if (nx > 0) {
                v[nx-1][j] = 0.0; // Tangencial en la pared derecha i=nx (afecta a la cara v más cercana)
            }
        }

        // --- NO-SLIP TANGENCIAL EN CARAS ADYACENTES A SÓLIDOS INTERNOS ---
        // Esto es más complejo de definir perfectamente solo con las caras u,v
        // sin afectar el flujo normal que ya fue manejado por enforceSolidCellBoundaryConditions().
        // Un enfoque es asegurar que las velocidades DENTRO de las celdas fluidas adyacentes a sólidos
        // tengan su componente tangencial a la pared sólida puesta a cero.
        // Sin embargo, esto es más naturalmente manejado en G2P o advección de partículas.

        // Por ahora, nos enfocaremos en el efecto a nivel de partícula en advectParticles,
        // ya que el no-slip a nivel de rejilla para sólidos internos puede ser complicado
        // y el `enforceSolidCellBoundaryConditions` ya maneja la componente normal.
        // Lo importante es que `useNoSlip` se respete en la colisión de partículas.
    }
    
    // Dentro de la clase FLIPSimulator:

    // Dentro de la clase FLIPSimulator, en el método projectPressure:
// Dentro de la clase FLIPSimulator:

// Dentro de la clase FLIPSimulator:
*/

applyTangentialNoSlip() {
    const grid = this.grid;
    const u = grid.u;
    const v = grid.v;
    const nx = grid.nx;
    const ny = grid.ny;

    if (!this.useNoSlip) {
        return;
    }

    // Paredes horizontales del dominio (superior e inferior) - afecta a 'u' (tangencial)
    for (let i = 0; i <= nx; i++) {
        if (!this.boundaryOpen.bottom) { // Solo si el contorno inferior está CERRADO
            u[i][0] = 0.0;
        }
        if (ny > 0 && !this.boundaryOpen.top) { // Solo si el contorno superior está CERRADO
             u[i][ny-1] = 0.0; 
        }
    }
    // Paredes verticales del dominio (izquierda y derecha) - afecta a 'v' (tangencial)
    for (let j = 0; j <= ny; j++) {
        if (!this.boundaryOpen.left) { // Solo si el contorno izquierdo está CERRADO
            v[0][j] = 0.0;
        }
        if (nx > 0 && !this.boundaryOpen.right) { // Solo si el contorno derecho está CERRADO
            v[nx-1][j] = 0.0; 
        }
    }
}

projectPressure(dt) {
    const grid = this.grid;
    const u = grid.u;
    const v = grid.v;
    const p = grid.p;
    const isFluid = grid.isFluid;
    const isSolid = grid.isSolid;
    const nx = grid.nx;
    const ny = grid.ny;
    const h = grid.cellSize;
    const simDensity = this.SIM_DENSITY;
    const relaxationFactor = this.pressureRelaxationFactor;
    const tolerance = this.pressureTolerance;
    const maxIterations = this.numPressureIterations;

    const targetDensity = this.targetCellDensity; // Usar el nuevo parámetro
    const densityCorrFactor = this.densityCorrectionFactor; // Usar el nuevo parámetro


    const divergence = Array(nx).fill(null).map(() => Array(ny).fill(0.0));
    const invH = 1.0 / h;

    // console.log(`En projectPressure - targetCellDensity: ${this.targetCellDensity.toFixed(2)}, densityCorrectionFactor: ${this.densityCorrectionFactor.toFixed(3)}, MAX_REASONABLE_CELL_DENSITY: ${this.MAX_REASONABLE_CELL_DENSITY.toFixed(2)}`);

    // --- 1. Calcular Divergencia ---
    for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
            if (!isFluid[i][j]) {
                divergence[i][j] = 0.0;
                continue;
            }
            // Divergencia estándar del campo de velocidades
            let currentDivergence = (u[i + 1][j] - u[i][j]) * invH +
                                      (v[i][j + 1] - v[i][j]) * invH;

            // NUEVO: Término de corrección por densidad
            const currentCellDensity = grid.cellDensity[i][j]; // Obtener la densidad de la rejilla
            currentDivergence -= densityCorrFactor * Math.max(0.0, (currentCellDensity - targetDensity));

            divergence[i][j] = currentDivergence;
        }
    }

    // --- 2. Resolver Presiones con Gauss-Seidel, Sub-Relajación y Manejo Correcto de BC ---
    const pressureRHSFactor = (simDensity * h * h) / dt;

    for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
            if (!isFluid[i][j] && !isSolid[i][j]) { 
                p[i][j] = 0.0;
            }
            // Si es FLUIDO, p[i][j] conserva su valor como conjetura inicial.
        }
    }

    // --- MODIFICACIÓN: Definir el límite para el clamp de presión ---
    const MAX_P_DELTA_PER_UPDATE = 1.0; // ¡VALOR DE PRUEBA! Ajústalo (0.05, 0.2, etc.)
    // Si es muy bajo, la presión no se propagará.
    // Si es muy alto, no tendrá efecto de clamp.
    // --- FIN MODIFICACIÓN ---


    let iter;
    for (iter = 0; iter < maxIterations; iter++) {
        let maxPressureChangeThisIteration = 0.0;

        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                if (!isFluid[i][j]) { 
                    continue; 
                }

                const p_old_for_cell = p[i][j]; 

                let sum_fluid_or_air_neighbor_pressures = 0;
                let num_solid_or_domain_boundary_faces = 0;

                // Vecino Derecho (P[i+1][j])
                if (i + 1 < nx) { 
                    if (isSolid[i + 1][j]) {
                        num_solid_or_domain_boundary_faces++;
                    } else { 
                        sum_fluid_or_air_neighbor_pressures += p[i + 1][j]; 
                    }
                } else { 
                    num_solid_or_domain_boundary_faces++;
                }

                // Vecino Izquierdo (P[i-1][j])
                if (i - 1 >= 0) { 
                    if (isSolid[i - 1][j]) {
                        num_solid_or_domain_boundary_faces++;
                    } else {
                        sum_fluid_or_air_neighbor_pressures += p[i - 1][j];
                    }
                } else { 
                    num_solid_or_domain_boundary_faces++;
                }

                // Vecino Superior (P[i][j+1])
                if (j + 1 < ny) { 
                    if (isSolid[i][j + 1]) {
                        num_solid_or_domain_boundary_faces++;
                    } else {
                        sum_fluid_or_air_neighbor_pressures += p[i][j + 1];
                    }
                } else { 
                    num_solid_or_domain_boundary_faces++;
                }

                // Vecino Inferior (P[i][j-1])
                if (j - 1 >= 0) { 
                    if (isSolid[i][j - 1]) {
                        num_solid_or_domain_boundary_faces++;
                    } else {
                        sum_fluid_or_air_neighbor_pressures += p[i][j - 1];
                    }
                } else { 
                    num_solid_or_domain_boundary_faces++;
                }
                
                const effective_diagonal_coefficient = 4.0 - num_solid_or_domain_boundary_faces;
                
                let p_candidate = p_old_for_cell; 
                if (effective_diagonal_coefficient > 0) {
                     p_candidate = (sum_fluid_or_air_neighbor_pressures - divergence[i][j] * pressureRHSFactor) 
                                   / effective_diagonal_coefficient;
                }
                
                // --- MODIFICACIÓN: Aplicar clamp de presión ---
                const p_after_relaxation = p_old_for_cell + relaxationFactor * (p_candidate - p_old_for_cell);
                const change_this_update = p_after_relaxation - p_old_for_cell;
                const clamped_change = Math.max(-MAX_P_DELTA_PER_UPDATE, Math.min(MAX_P_DELTA_PER_UPDATE, change_this_update));
                const p_new_for_cell = p_old_for_cell + clamped_change;
                // --- FIN MODIFICACIÓN ---
                
                maxPressureChangeThisIteration = Math.max(maxPressureChangeThisIteration, Math.abs(clamped_change)); // Usar el cambio real
                p[i][j] = p_new_for_cell;
            }
        }

        if (maxPressureChangeThisIteration < tolerance && iter > 0) {
            iter++; break; 
        }
    }
    const actualIterationsPerformed = (maxIterations === 0) ? 0 : iter;

    // --- 3. Actualizar Velocidades ---
    const velocityCorrectionFactor = dt / (simDensity * h);

    // Actualizar velocidades U (caras verticales)
    // u[i][j] es la cara entre P[i-1][j] y P[i][j]
    for (let i = 1; i < nx; i++) { // Itera sobre las caras u INTERNAS del dominio
        for (let j = 0; j < ny; j++) {
            // Solo actualizar si esta cara NO es una frontera con un sólido
            // y al menos una de las celdas adyacentes es fluido.
            if (!isSolid[i - 1][j] && !isSolid[i][j]) { // La cara está entre dos celdas NO SÓLIDAS
                if (isFluid[i - 1][j] || isFluid[i][j]) { // Y al menos una es fluido (o ambas son aire)
                    u[i][j] -= velocityCorrectionFactor * (p[i][j] - p[i - 1][j]);
                }
            }
            // Si una de las celdas es sólida (ej: isSolid[i-1][j]), entonces u[i][j] es una
            // velocidad de contorno que ya fue puesta a 0 y no debe cambiar.
            // La condición de Neumann en el solver de presión debería hacer que p[i][j] - p[i-1][j] (fantasma) sea 0.
        }
    }

    // Actualizar velocidades V (caras horizontales)
    // v[i][j] es la cara entre P[i][j-1] y P[i][j]
    for (let i = 0; i < nx; i++) {
        for (let j = 1; j < ny; j++) { // Itera sobre las caras v INTERNAS del dominio
            // Solo actualizar si esta cara NO es una frontera con un sólido
            // y al menos una de las celdas adyacentes es fluido.
            if (!isSolid[i][j - 1] && !isSolid[i][j]) { // La cara está entre dos celdas NO SÓLIDAS
                if (isFluid[i][j - 1] || isFluid[i][j]) { // Y al menos una es fluido (o ambas son aire)
                    v[i][j] -= velocityCorrectionFactor * (p[i][j] - p[i][j - 1]);
                }
            }
            // Si una de las celdas es sólida (ej: isSolid[i][j-1]), entonces v[i][j] es una
            // velocidad de contorno que ya fue puesta a 0 y no debe cambiar.
        }
    }
    return actualIterationsPerformed;
} // Fin de projectPressure

    _interpolateFromUGrid(px, py, uGrid) {
        // ... (Código completo de _interpolateFromUGrid como lo teníamos) ...
        const h = this.grid.cellSize; const nx = this.grid.nx; const ny = this.grid.ny; const gx = px/h; const gy = py/h - 0.5; const i0 = Math.floor(gx); const j0 = Math.floor(gy); const tx = gx - i0; const ty = gy - j0;
        const u00 = (i0 >= 0 && i0 <= nx && j0 >= 0 && j0 < ny) ? uGrid[i0][j0] : 0; const u10 = (i0+1 >= 0 && i0+1 <= nx && j0 >= 0 && j0 < ny) ? uGrid[i0+1][j0] : 0; const u01 = (i0 >= 0 && i0 <= nx && j0+1 >= 0 && j0+1 < ny) ? uGrid[i0][j0+1] : 0; const u11 = (i0+1 >= 0 && i0+1 <= nx && j0+1 >= 0 && j0+1 < ny) ? uGrid[i0+1][j0+1] : 0;
        return (1-tx)*(1-ty)*u00 + tx*(1-ty)*u10 + (1-tx)*ty*u01 + tx*ty*u11;
    }

    _interpolateFromVGrid(px, py, vGrid) {
        // ... (Código completo de _interpolateFromVGrid como lo teníamos) ...
        const h = this.grid.cellSize; const nx = this.grid.nx; const ny = this.grid.ny; const gx = px/h - 0.5; const gy = py/h; const i0 = Math.floor(gx); const j0 = Math.floor(gy); const tx = gx - i0; const ty = gy - j0;
        const v00 = (i0 >= 0 && i0 < nx && j0 >= 0 && j0 <= ny) ? vGrid[i0][j0] : 0; const v10 = (i0+1 >= 0 && i0+1 < nx && j0 >= 0 && j0 <= ny) ? vGrid[i0+1][j0] : 0; const v01 = (i0 >= 0 && i0 < nx && j0+1 >= 0 && j0+1 <= ny) ? vGrid[i0][j0+1] : 0; const v11 = (i0+1 >= 0 && i0+1 < nx && j0+1 >= 0 && j0+1 <= ny) ? vGrid[i0+1][j0+1] : 0;
        return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11;
    }

    gridToParticle(alpha) {
        // ... (Código completo de gridToParticle como lo teníamos) ...
        const gridFinal_u = this.grid.u; const gridFinal_v = this.grid.v; const gridBefore_u = this.grid.u_before_pressure; const gridBefore_v = this.grid.v_before_pressure;
        for (const particle of this.particles) { const oldParticleVx = particle.vx; const oldParticleVy = particle.vy; const u_pic = this._interpolateFromUGrid(particle.x, particle.y, gridFinal_u); const v_pic = this._interpolateFromVGrid(particle.x, particle.y, gridFinal_v); const u_grid_before = this._interpolateFromUGrid(particle.x, particle.y, gridBefore_u); const v_grid_before = this._interpolateFromVGrid(particle.x, particle.y, gridBefore_v); const delta_u_grid = u_pic - u_grid_before; const delta_v_grid = v_pic - v_grid_before; const flipVx = oldParticleVx + delta_u_grid; const flipVy = oldParticleVy + delta_v_grid; particle.vx = alpha * flipVx + (1-alpha) * u_pic; particle.vy = alpha * flipVy + (1-alpha) * v_pic; }
    }


    advectParticles(dt) {
        const particles = this.particles;
        const h = this.grid.cellSize; // No se usa directamente aquí si tenemos AABBs, pero epsilon sí
        const grid = this.grid; // Necesario para isSolid, nx, ny, cellSize
        const nx = this.grid.nx; // Necesario para los límites del dominio
        const ny = this.grid.ny; // Necesario para los límites del dominio

        const domainBounds = {
            xMin: 0.0, xMax: nx * h,
            yMin: 0.0, yMax: ny * h
        };

        const restitution = this.restitutionCoefficient;
        const friction = this.frictionFactor;
        const epsilonClosedBoundary = 1e-6 * h; // Epsilon solo para contornos CERRADOS

        // --- SUB-STEPPING ---
        let numSubsteps = 2; // Por defecto, 1 sub-paso (sin sub-stepping real)
        if (this.hasComplexSolids) {
            numSubsteps = 3; // O el número de sub-pasos que necesites para diagonales
                            // Podrías hacerlo configurable: this.substepsForComplex = 3;
        }
        const subDt = dt / numSubsteps;

            // Iterar sobre las partículas. Es más seguro marcar para eliminar y filtrar después.
        for (const particle of this.particles) {
            particle.toBeRemoved = false; // Inicializar/resetear el flag para este paso de advección

            for (let step = 0; step < numSubsteps; step++) {
                if (particle.toBeRemoved) break; // Si ya está marcada, no seguir con subpasos para esta partícula

                const oldX_substep = particle.x;
                const oldY_substep = particle.y;
                const originalVx_substep = particle.vx;
                const originalVy_substep = particle.vy;

                // --- Advección y Colisión en Eje X ---
                let tentativeX = particle.x + originalVx_substep * subDt;

                if (tentativeX <= domainBounds.xMin) {
                    if (this.boundaryOpen.left) {
                        particle.toBeRemoved = true;
                    } else { // Contorno izquierdo cerrado
                        tentativeX = domainBounds.xMin + epsilonClosedBoundary;
                        particle.vx = -restitution * originalVx_substep;
                        if (this.useNoSlip) particle.vy = 0.0; else particle.vy *= friction;
                    }
                } else if (tentativeX >= domainBounds.xMax) {
                    if (this.boundaryOpen.right) {
                        particle.toBeRemoved = true;
                    } else { // Contorno derecho cerrado
                        tentativeX = domainBounds.xMax - epsilonClosedBoundary;
                        particle.vx = -restitution * originalVx_substep;
                        if (this.useNoSlip) particle.vy = 0.0; else particle.vy *= friction;
                    }
                }
                particle.x = tentativeX; // Actualizar posición X tentativa o final del subpaso
                if (particle.toBeRemoved) continue; // Si se marcó para eliminar, saltar resto del subpaso

                // Colisión con sólidos internos en X (esta lógica se mantiene)
                const cjForXCheck = Math.floor(oldY_substep / h);
                const ciNewX = Math.floor(particle.x / h); // Usar particle.x ya actualizado
                if (ciNewX >= 0 && ciNewX < nx && cjForXCheck >= 0 && cjForXCheck < ny &&
                    grid.isSolid[ciNewX][cjForXCheck]) {
                    const ciOldX = Math.floor(oldX_substep / h);
                    if (!grid.isSolid[ciOldX][cjForXCheck] || ciOldX !== ciNewX) {
                        if (originalVx_substep > 0) {
                            particle.x = ciNewX * h - epsilonClosedBoundary;
                        } else if (originalVx_substep < 0) {
                            particle.x = (ciNewX + 1) * h + epsilonClosedBoundary;
                        }
                        // Aplicar respuesta de velocidad solo si realmente hubo una colisión que revirtió el movimiento
                        if ((originalVx_substep > 0 && particle.x < oldX_substep + originalVx_substep * subDt) ||
                            (originalVx_substep < 0 && particle.x > oldX_substep + originalVx_substep * subDt) ||
                            originalVx_substep === 0) {
                            particle.vx = -restitution * originalVx_substep;
                            if (this.useNoSlip) { particle.vy = 0.0; } 
                            else { particle.vy *= friction; }
                        }
                    }
                }

                // --- Advección y Colisión en Eje Y ---
                const vyForYSubstep = particle.vy; // Usar vy que pudo ser afectada por colisión X
                let tentativeY = oldY_substep + vyForYSubstep * subDt;

                if (tentativeY <= domainBounds.yMin) {
                    if (this.boundaryOpen.bottom) {
                        particle.toBeRemoved = true;
                    } else { // Contorno inferior cerrado
                        tentativeY = domainBounds.yMin + epsilonClosedBoundary;
                        particle.vy = -restitution * vyForYSubstep;
                        if (this.useNoSlip) particle.vx = 0.0; else particle.vx *= friction;
                    }
                } else if (tentativeY >= domainBounds.yMax) {
                    if (this.boundaryOpen.top) {
                        particle.toBeRemoved = true;
                    } else { // Contorno superior cerrado
                        tentativeY = domainBounds.yMax - epsilonClosedBoundary;
                        particle.vy = -restitution * vyForYSubstep;
                        if (this.useNoSlip) particle.vx = 0.0; else particle.vx *= friction;
                    }
                }
                particle.y = tentativeY; // Actualizar posición Y tentativa o final del subpaso
                if (particle.toBeRemoved) continue; // Si se marcó para eliminar, saltar resto del subpaso

                // Colisión con sólidos internos en Y (esta lógica se mantiene)
                const ciForYCheck = Math.floor(particle.x / h); // Usar particle.x ya actualizado
                const cjNewY = Math.floor(particle.y / h); // Usar particle.y ya actualizado
                if (ciForYCheck >= 0 && ciForYCheck < nx && cjNewY >= 0 && cjNewY < ny &&
                    grid.isSolid[ciForYCheck][cjNewY]) {
                    const cjOldY = Math.floor(oldY_substep / h);
                    if (!grid.isSolid[ciForYCheck][cjOldY] || cjNewY !== cjOldY) {
                        if (vyForYSubstep > 0) {
                            particle.y = cjNewY * h - epsilonClosedBoundary;
                        } else if (vyForYSubstep < 0) {
                            particle.y = (cjNewY + 1) * h + epsilonClosedBoundary;
                        }
                        if ((vyForYSubstep > 0 && particle.y < oldY_substep + vyForYSubstep * subDt) ||
                            (vyForYSubstep < 0 && particle.y > oldY_substep + vyForYSubstep * subDt) ||
                            vyForYSubstep === 0) {
                            particle.vy = -restitution * vyForYSubstep;
                            if (this.useNoSlip) { particle.vx = 0.0; } 
                            else { particle.vx *= friction; }
                        }
                    }
                }
            } // Fin del bucle de sub-pasos
        } // Fin del bucle de partículas

        // Filtrar las partículas marcadas para eliminación DESPUÉS de iterar por todas
        const oldParticleCount = this.particles.length;
        this.particles = this.particles.filter(p => !p.toBeRemoved);
        const removedCount = oldParticleCount - this.particles.length;
        // if (removedCount > 0) {
        //     console.log(`Eliminadas ${removedCount} partículas por contornos abiertos.`);
        // }
        // / console.log("Módulo: Advección de Partículas Completado (con sub-pasos).");
    }

    step(dt) {
        this.simulationTime += dt; 

        this.particleToGrid();
        for (let i = 0; i <= this.grid.nx; i++) { for (let j = 0; j < this.grid.ny; j++) { this.grid.u_before_pressure[i][j] = this.grid.u[i][j]; }}
        for (let i = 0; i < this.grid.nx; i++) { for (let j = 0; j <= this.grid.ny; j++) { this.grid.v_before_pressure[i][j] = this.grid.v[i][j]; }}
        
        this.applyExternalForces(dt);
        this.enforceMACBoundaryConditions(); 
        
        this.enforceSolidCellBoundaryConditions(); // Condiciones de contorno para SÓLIDOS INTERNOS


        const actualPressureIterations = this.projectPressure(dt); // <<--- ALMACENAR RESULTADO

        //this.projectPressure(dt);
        
        /*if (this.useNoSlip) {
            this.applyTangentialNoSlip();
        }*/
        
        this.applyTangentialNoSlip(); // APLICA A LÍMITES DEL DOMINIO (principalmente)
                                // El no-slip en sólidos internos se maneja en advectParticles

        this.gridToParticle(this.picFlipAlpha);
        this.advectParticles(dt);

        if (this.particles.length > 0) {
            let currentMaxHeight = 0.0;
            for (const particle of this.particles) {
                if (particle.y > currentMaxHeight) {
                    currentMaxHeight = particle.y;
                }
            }
            this.maxHeightData.push({ time: this.simulationTime, height: currentMaxHeight });
            //const startTimeWindow = this.simulationTime - this.graphTimeWindowSeconds;
            //this.maxHeightData = this.maxHeightData.filter(point => point.time >= startTimeWindow);
        }
        // --- NUEVO: Registrar número de iteraciones de presión ---
        this.pressureIterationsData.push({ time: this.simulationTime, iterations: actualPressureIterations });
        // --- FIN NUEVO ---
        this.particleCountHistoryData.push({ time: this.simulationTime, count: this.particles.length });

    }

    // ESTE ES EL NUEVO MÉTODO QUE SE AÑADE
    getMaxHeightHistory() {
        return this.maxHeightData;
    }

    getPressureIterationsHistory() {
        return this.pressureIterationsData;
    }

    getParticleCountHistory() {
        return this.particleCountHistoryData;
    }
} // Fin de la clase FLIPSimulator

// --- Ejemplo de cómo podrías empezar a usarlo (sin cambios) ---
// const simWidth = 2.0; 
// const simHeight = 1.0; 
// const cellResolution = 0.05; 
// const simulator = new FLIPSimulator(simWidth, simHeight, cellResolution);
// simulator.initializeFluidVolume(0.0, 0.0, simWidth * 0.5, simHeight * 0.75, 2);
// const deltaTime = 0.016; 
// simulator.step(deltaTime); // Esto ahora llamará a tu particleToGrid() implementado