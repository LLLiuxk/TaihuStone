import * as THREE from 'three';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// ============================================================
// Trial Configuration
// ============================================================
const trials = [
    {
        id: 1,
        question: "Which structure feels more 'open' and 'breathable'?",
        type: "VP Evaluation",
        modelA: "models/dummyA.stl",
        modelB: "models/dummyB.stl"
    },
    {
        id: 2,
        question: "Which structure looks more like a naturally weathered stone?",
        type: "Organic Likeness",
        modelA: "models/dummyE.stl",
        modelB: "models/dummyF.stl"
    }
    // Add more trials here...
];

// ============================================================
// State
// ============================================================
let currentTrialIndex = 0;
const results = [];
let isAutoRotating = true;
let isSyncRotation = false;
let _syncLock = false; // prevents infinite loop when syncing

// Initial camera position — change this to adjust default view angle
const INIT_CAM = { x: 2, y: 0.5, z: 2 };

// We keep references to the two viewer instances so we can dispose/rebuild
let viewerInstanceA = null;
let viewerInstanceB = null;

// ============================================================
// DOM references
// ============================================================
const qText    = document.getElementById('question-text');
const pText    = document.getElementById('progress-text');
const pFill    = document.getElementById('progress-fill');
const containerA = document.getElementById('viewer-a');
const containerB = document.getElementById('viewer-b');
const trialMain  = document.getElementById('trial-container');
const compMain   = document.getElementById('completion-container');
const resJson    = document.getElementById('results-json');
const btnPause   = document.getElementById('btn-pause');
const btnSync    = document.getElementById('btn-sync');

document.getElementById('btn-a').addEventListener('click', () => selectModel('A'));
document.getElementById('btn-b').addEventListener('click', () => selectModel('B'));
document.getElementById('btn-download').addEventListener('click', downloadResults);
document.getElementById('btn-reset').addEventListener('click', resetViews);
btnPause.addEventListener('click', toggleRotation);
btnSync.addEventListener('click', toggleSync);

// ============================================================
// Control functions
// ============================================================
function resetViews() {
    [viewerInstanceA, viewerInstanceB].forEach(v => {
        if (!v) return;
        v.camera.position.set(INIT_CAM.x, INIT_CAM.y, INIT_CAM.z);
        v.controls.target.set(0, 0, 0);
        v.controls.update();
    });
}

function toggleRotation() {
    isAutoRotating = !isAutoRotating;
    [viewerInstanceA, viewerInstanceB].forEach(v => {
        if (!v) return;
        v.controls.autoRotate = isAutoRotating;
    });
    if (isAutoRotating) {
        btnPause.innerHTML = '<span class="ctrl-icon">\u23f8</span> Pause Rotation';
        btnPause.classList.remove('active');
    } else {
        btnPause.innerHTML = '<span class="ctrl-icon">\u25b6</span> Resume Rotation';
        btnPause.classList.add('active');
    }
}

function toggleSync() {
    isSyncRotation = !isSyncRotation;
    if (isSyncRotation) {
        btnSync.innerHTML = '<span class="ctrl-icon">\ud83d\udd17</span> Synced';
        btnSync.classList.add('active');
        // Wire up listeners
        setupSyncListeners();
        // Prevent double auto-rotation speed by disabling it on viewer B
        if (viewerInstanceB) {
            viewerInstanceB.controls.autoRotate = false;
        }
    } else {
        btnSync.innerHTML = '<span class="ctrl-icon">\ud83d\udd17</span> Sync Rotation';
        btnSync.classList.remove('active');
        // Remove listeners
        removeSyncListeners();
        // Restore auto-rotation state
        if (viewerInstanceB) {
            viewerInstanceB.controls.autoRotate = isAutoRotating;
        }
    }
}

// --- Sync helpers ---
function onControlsChangeA() {
    if (_syncLock || !isSyncRotation || !viewerInstanceB) return;
    _syncLock = true;
    mirrorCamera(viewerInstanceA, viewerInstanceB);
    _syncLock = false;
}

function onControlsChangeB() {
    if (_syncLock || !isSyncRotation || !viewerInstanceA) return;
    _syncLock = true;
    mirrorCamera(viewerInstanceB, viewerInstanceA);
    _syncLock = false;
}

function mirrorCamera(src, dst) {
    dst.camera.position.copy(src.camera.position);
    dst.controls.target.copy(src.controls.target);
    dst.controls.update();
}

function setupSyncListeners() {
    if (viewerInstanceA) viewerInstanceA.controls.addEventListener('change', onControlsChangeA);
    if (viewerInstanceB) viewerInstanceB.controls.addEventListener('change', onControlsChangeB);
}

function removeSyncListeners() {
    if (viewerInstanceA) viewerInstanceA.controls.removeEventListener('change', onControlsChangeA);
    if (viewerInstanceB) viewerInstanceB.controls.removeEventListener('change', onControlsChangeB);
}

// ============================================================
// Three.js viewer factory
// ============================================================
function createViewer(container) {
    // Clean up any previous content
    while (container.firstChild) container.removeChild(container.firstChild);

    const width  = container.clientWidth  || 500;
    const height = container.clientHeight || 500;

    // Scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xf9f9f9);

    // Camera
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 1000);
    camera.position.set(2, 2, 2);

    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.shadowMap.enabled = true;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.2;
    container.appendChild(renderer.domElement);

    // Controls (interactive drag-rotate + auto-rotate)
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = false; // Disabling damping makes pausing and syncing instantaneous
    controls.autoRotate = true;
    controls.autoRotateSpeed = 2.0;

    // Lighting — key + fill + rim + ambient
    const keyLight = new THREE.DirectionalLight(0xffffff, 2.5);
    keyLight.position.set(5, 8, 5);
    keyLight.castShadow = true;
    scene.add(keyLight);

    const fillLight = new THREE.DirectionalLight(0x8ec8f0, 1.0);
    fillLight.position.set(-4, 3, -3);
    scene.add(fillLight);

    const rimLight = new THREE.DirectionalLight(0xfff4e6, 0.8);
    rimLight.position.set(0, -2, -5);
    scene.add(rimLight);

    const ambient = new THREE.AmbientLight(0x404060, 0.6);
    scene.add(ambient);

    // Subtle ground plane for shadow / depth cue
    const planeGeo = new THREE.PlaneGeometry(10, 10);
    const planeMat = new THREE.ShadowMaterial({ opacity: 0.25 });
    const plane = new THREE.Mesh(planeGeo, planeMat);
    plane.rotation.x = -Math.PI / 2;
    plane.position.y = -1;
    plane.receiveShadow = true;
    scene.add(plane);

    // Animation loop
    let animId;
    function animate() {
        animId = requestAnimationFrame(animate);
        controls.update();
        renderer.render(scene, camera);
    }
    animate();

    // Resize handling
    const onResize = () => {
        const w = container.clientWidth;
        const h = container.clientHeight;
        camera.aspect = w / h;
        camera.updateProjectionMatrix();
        renderer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    return {
        scene, camera, controls, renderer,
        dispose() {
            cancelAnimationFrame(animId);
            window.removeEventListener('resize', onResize);
            controls.dispose();
            renderer.dispose();
            while (container.firstChild) container.removeChild(container.firstChild);
        }
    };
}

// ============================================================
// Load an STL into an existing viewer
// ============================================================
function loadSTL(viewer, url) {
    return new Promise((resolve, reject) => {
        const loader = new STLLoader();
        loader.load(
            url,
            (geometry) => {
                // Convert from Z-up (CAD/3D-printing convention) to Y-up (Three.js convention)
                geometry.rotateX(-Math.PI / 2);

                // Centre the geometry
                geometry.computeBoundingBox();
                const box = geometry.boundingBox;
                const center = new THREE.Vector3();
                box.getCenter(center);
                geometry.translate(-center.x, -center.y, -center.z);

                // Scale to fit in a ~2-unit cube
                const size = new THREE.Vector3();
                box.getSize(size);
                const maxDim = Math.max(size.x, size.y, size.z);
                const scale = 2.0 / maxDim;
                geometry.scale(scale, scale, scale);

                // Compute normals for smooth shading
                geometry.computeVertexNormals();

                // Clay-like material
                const material = new THREE.MeshStandardMaterial({
                    color: 0xd4c5a9,
                    roughness: 0.65,
                    metalness: 0.05,
                    flatShading: false
                });

                const mesh = new THREE.Mesh(geometry, material);
                mesh.castShadow = true;
                mesh.receiveShadow = true;
                viewer.scene.add(mesh);

                // Adjust camera to frame the model
                viewer.camera.position.set(INIT_CAM.x, INIT_CAM.y, INIT_CAM.z);
                viewer.controls.target.set(0, 0, 0);
                viewer.controls.autoRotate = isAutoRotating;
                viewer.controls.update();

                resolve(mesh);
            },
            undefined,
            (err) => {
                console.error('STL load error:', err);
                reject(err);
            }
        );
    });
}

// ============================================================
// Trial management
// ============================================================
function initTrial() {
    if (currentTrialIndex >= trials.length) {
        showCompletion();
        return;
    }

    const trial = trials[currentTrialIndex];

    // Update UI text
    qText.textContent = trial.question;
    pText.textContent = `Trial ${currentTrialIndex + 1} / ${trials.length}`;
    pFill.style.width = `${(currentTrialIndex / trials.length) * 100}%`;

    // Dispose old viewers
    if (viewerInstanceA) viewerInstanceA.dispose();
    if (viewerInstanceB) viewerInstanceB.dispose();

    // Create fresh viewers
    viewerInstanceA = createViewer(containerA);
    viewerInstanceB = createViewer(containerB);

    // Load models
    loadSTL(viewerInstanceA, trial.modelA);
    loadSTL(viewerInstanceB, trial.modelB);
}

function selectModel(choice) {
    const trial = trials[currentTrialIndex];
    results.push({
        trialId: trial.id,
        questionType: trial.type,
        modelA: trial.modelA,
        modelB: trial.modelB,
        selected: choice === 'A' ? trial.modelA : trial.modelB,
        choiceLabel: choice,
        timestamp: new Date().toISOString()
    });
    currentTrialIndex++;
    initTrial();
}

function showCompletion() {
    pFill.style.width = '100%';
    pText.textContent = 'Completed';
    trialMain.style.display = 'none';
    compMain.style.display = 'block';
    resJson.value = JSON.stringify(results, null, 2);
}

function downloadResults() {
    const blob = new Blob([JSON.stringify(results, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `taihustone_results_${Date.now()}.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

// ============================================================
// Boot
// ============================================================
initTrial();
