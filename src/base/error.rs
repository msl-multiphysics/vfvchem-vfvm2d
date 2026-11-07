use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum VFVChemError {
    // vertices
    #[error("Points: Maximum sampling attempts `max_try` must be greater than zero. Got max_try = {max_try}.")]
    PointsNonPositiveAttempts {
        max_try: usize,
    },
    #[error("Points: Spacing `h` must be finite and at least equal to the node spacing. Got h = {h}.")]
    PointsInvalidSpacing {
        h: f64,
    },
    #[error("Points: Ramp length `ramp` must be positive and finite. Got ramp = {ramp}.")]
    PointsNonPositiveRampLength {
        ramp: f64,
    },
    #[error("Nodes: x_max must be greater than x_min. Got x_min = {x_min} and x_max = {x_max}.")]
    NodesInvalidBoundsX {
        x_min: f64,
        x_max: f64,
    },
    #[error("Nodes: y_max must be greater than y_min. Got y_min = {y_min} and y_max = {y_max}.")]
    NodesInvalidBoundsY {
        y_min: f64,
        y_max: f64,
    },
    #[error("Nodes: Spacing `h` must be positive and finite. Got h = {h}.")]
    NodesNonPositiveSpacing {
        h: f64,
    },

    // geometry
    #[error("Domain: Group ID {grp2d_id} not found.")]
    DomainGroupIDNotFound {
        grp2d_id: String,
    },
    #[error("Boundary: Group ID {grp1d_id} not found in Nodes.")]
    BoundaryGroupIDNotFound {
        grp1d_id: String,
    },
    #[error("Interface: Group ID {grp1d_id} not found in Nodes.")]
    InterfaceGroupIDNotFound {
        grp1d_id: String,
    },
    #[error("Interface: Domains {dom_a} and {dom_b} do not share an interface at group ID {grp1d_id}.")]
    InterfaceNoCommonInterface {
        dom_a: usize,
        dom_b: usize,
        grp1d_id: String,
    },
    #[error("Triangle: Error encountered in Delaunay triangulation: {source}")]
    TriangleError {
        #[source]
        source: anyhow::Error,
    },

    // variables
    #[error("ScalarDomain: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    ScalarDomainFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },
    #[error("ScalarBoundary: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    ScalarBoundaryFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },
    #[error("ScalarInterface: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    ScalarInterfaceFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },
    #[error("VectorDomain: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    VectorDomainFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },
    #[error("VectorBoundary: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    VectorBoundaryFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },
    #[error("VectorInterface: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    VectorInterfaceFunctionDependsOnNonUnknown {
        scldom_id: usize,
    },

    // read and write
    #[error("ReadObj: Could not read file '{path}': {source}")]
    ReadObjError {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("ReadObj: Error encountered while parsing file '{path}': {message}")]
    ReadObjParse {
        path: String,
        message: String,
    },
    #[error("WriteVtu: Could not write file '{path}': {source}")]
    WriteVtuError {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("WriteCsv: Could not write file '{path}': {source}")]
    WriteCsvError {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },

    // linear
    #[error("Linear: Number of threads `num_thread` must be positive. Got num_thread = {num_thread}")]
    LinearNonPositiveThreads {
        num_thread: usize,
    },
    #[error("Linear: Error encountered in rslab solver: {source}")]
    LinearRslabError {
        #[source]
        source: rslab::RslabError,
    },

    // operators
    #[error("Operator: Use apply_time instead of apply.")]
    OperatorUseApplyTime {},
    #[error("Operator: Only unknown-type scalars and vectors can be solved for in the global linear system.")]
    OperatorNonUnknownInLinearSystem {},
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarDomainNotOnDomain {
        scldom_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarDomainNotOnBoundary {
        scldom_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarBoundary ID {sclbnd_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarBoundaryNotOnBoundary {
        sclbnd_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarDomain IDs {scldom_a_id} and {scldom_b_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarDomainNotOnInterface {
        scldom_a_id: usize,
        scldom_b_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarInterface ID {sclitf_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarInterfaceNotOnInterface {
        sclitf_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarFlowNotOnDomain {
        sclflw_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarFlowNotOnBoundary {
        sclflw_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorDomain ID {vecdom_id} is not defined on Domain ID {dom_id}.")]
    OperatorVectorDomainNotOnDomain {
        vecdom_id: usize,
        dom_id: usize,
    },
    #[error("Operator: VectorDomain ID {vecdom_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorVectorDomainNotOnBoundary {
        vecdom_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorBoundary ID {vecbnd_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorVectorBoundaryNotOnBoundary {
        vecbnd_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorDomain IDs {vecdom_a_id} and {vecdom_b_id} is not defined on Interface ID {itf_id}.")]
    OperatorVectorDomainNotOnInterface {
        vecdom_a_id: usize,
        vecdom_b_id: usize,
        itf_id: usize,
    },
    #[error("Operator: VectorInterface ID {vecitf_id} is not defined on Interface ID {itf_id}.")]
    OperatorVectorInterfaceNotOnInterface {
        vecitf_id: usize,
        itf_id: usize,
    },

    // solvers
    #[error("Solver: Maximum iterations `max_iter` must be greater than zero. Got max_iter = {max_iter}.")]
    SolverNonPositiveMaxIter {
        max_iter: usize,
    },
    #[error("Solver: Tolerance `tol` must be greater than zero. Got tol = {tol}.")]
    SolverNonPositiveTolerance {
        tol: f64,
    },
    #[error("Solver: Damping `damp` must be in (0, 1]. Got damp = {damp}.")]
    SolverInvalidDamping {
        damp: f64,
    },
    #[error("Solver: Nonlinear iteration diverged.")]
    SolverNonFiniteResidual {},
    #[error("Solver: Nonlinear iteration did not converge within {max_iter} iterations.")]
    SolverDidNotConverge {
        max_iter: usize,
    },
    #[error("Solver: Time step `dt` must be greater than zero. Got dt = {dt}.")]
    SolverNonPositiveTimeStep {
        dt: f64,
    },
    #[error("Solver: Number of time steps `num_ts` must be greater than zero. Got num_ts = {num_ts}.")]
    SolverNonPositiveNumTimeSteps {
        num_ts: usize,
    },
    #[error("Solver: Write frequency `num_ts_write` must be greater than zero. Got num_ts_write = {num_ts_write}.")]
    SolverNonPositiveWriteFrequency {
        num_ts_write: usize,
    },
    #[error("Solver: Initial time step `dt_init` must be greater than zero. Got dt_init = {dt_init}.")]
    SolverNonPositiveInitialStepSize {
        dt_init: f64,
    },
    #[error("Solver: Minimum time step `dt_min` must be greater than the initial timestep. Got dt_min = {dt_min}.")]
    SolverInvalidMinimumStepSize {
        dt_min: f64,
    },
    #[error("Solver: Maximum time step `dt_max` must be greater than the minimum timestep. Got dt_max = {dt_max}.")]
    SolverInvalidMaximumStepSize {
        dt_max: f64,
    },
    #[error("Solver: Minimum growth factor `grow_min` must be in (0, 1]. Got grow_min = {grow_min}.")]
    SolverInvalidGrowthMinimum {
        grow_min: f64,
    },
    #[error("Solver: Maximum growth factor `grow_max` must be greater than one. Got grow_max = {grow_max}.")]
    SolverInvalidGrowthMaximum {
        grow_max: f64,
    },
    #[error("Solver: Growth factor rate `grow_rate` must be greater than one. Got grow_rate = {grow_rate}.")]
    SolverInvalidGrowthRate {
        grow_rate: f64,
    },
    #[error("Solver: Adaptive solver did not converge at minimum timestep {dt_min}.")]
    SolverAdaptiveDidNotConverge {
        dt_min: f64,
    },
    #[error("Solver: Output times must be non-negative and monotonically increasing, and the first step must be greater than dt_init.")]
    SolverAdaptiveInvalidOutputTime {},
}
