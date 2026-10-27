use vfvchem_vfvm2d::*;
use std::fs::{read_to_string, remove_dir_all};
use std::path::PathBuf;

// transient heat on a unit square: T = 0 initially, T = 0 on the left and right, insulated top and bottom
// heat source Q switched on for t in (t0, t0 + w]; constant Q if pulse is None
struct Heat { vars: Variables, phys: TransientHeat, temp: usize }

fn heat(pulse: Option<(f64, f64)>) -> Result<Heat, VFVChemError> {
    let mut vars = Variables::new_bounds(0.0, 0.0, 1.0, 1.0, 0.1, 0.1, 0.1, 100, 42)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;
    let temp = vars.add_scldom_unkn(dom, 0.0, "")?;
    let vlcp = vars.add_scldom_cnst(dom, 1.0, "")?;
    let cond = vars.add_scldom_cnst(dom, 0.1, "")?;
    let hsrc_func: SclDomFunc = match pulse {
        Some((t0, w)) => Box::new(move |t, _p, _u| if t > t0 && t <= t0 + w + 1e-12 { 100.0 } else { 0.0 }),
        None => Box::new(|_t, _p, _u| 10.0),
    };
    let hprd = vars.add_scldom_cnst(dom, 0.0, "")?;
    let hsrc = vars.add_scldom_func(dom, hsrc_func, vec![], "")?;
    let temp_w = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let temp_e = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let hflx_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;
    let hflx_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hprd, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_w)?;
    phys.add_temp_bnd(bnd_r, temp_e)?;
    phys.add_hflx_bnd(bnd_b, hflx_b)?;
    phys.add_hflx_bnd(bnd_t, hflx_t)?;
    Ok(Heat { vars, phys, temp })
}

fn folder(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("vfvchem_adaptive_{name}_{}", std::process::id()));
    let _ = remove_dir_all(&path);
    path
}

fn read_csv(path: PathBuf) -> (String, Vec<Vec<f64>>) {
    let text = read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap().to_string();
    let rows = lines.map(|l| l.split(',').map(|x| x.parse::<f64>().unwrap()).collect()).collect();
    (header, rows)
}

#[test]
fn options_defaults() {
    let opt = AdaptiveOptions::new(0.01, 2.0, "out").unwrap();
    assert_eq!(opt.dt_init, 0.01);
    assert_eq!(opt.t_end, 2.0);
    assert_eq!(opt.rtol, 1e-3);
    assert_eq!(opt.folder_path, "out");
    assert_eq!(opt.out_times.len(), 11);
    for (i, &t) in opt.out_times.iter().enumerate() {
        assert!((t - 0.2 * i as f64).abs() < 1e-12);
    }
}

#[test]
fn default_output_times() {
    let path = folder("default");
    let mut h = heat(None).unwrap();
    let mut lin: Box<dyn LinearBase> = FaerLu::new(1).unwrap();
    let mut opt = AdaptiveOptions::new(0.01, 1.0, path.to_str().unwrap()).unwrap();
    opt.atol = 1e-6;  // T starts at 0; a tiny atol would drive the first step to dt_min
    h.phys.solve_adaptive(&mut h.vars, &mut lin, &opt).unwrap();

    // file index i is written at t = 0.1 * i
    let (header, rows) = read_csv(path.join("adaptive_times.csv"));
    assert_eq!(header, "file_index,time");
    assert_eq!(rows.len(), 11);
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row[0], i as f64);
        assert!((row[1] - 0.1 * i as f64).abs() < 1e-12);
    }

    // one row per accepted step; summary counts repeated; steps end exactly at t_end
    let (header, rows) = read_csv(path.join("adaptive_stats.csv"));
    assert_eq!(header, "accepted,rejected_error,rejected_picard,picard_iterations,time,dt,error");
    assert_eq!(rows.len(), rows[0][0] as usize);
    assert!(rows.iter().all(|r| r[..4] == rows[0][..4]));
    let dt_sum: f64 = rows.iter().map(|r| r[5]).sum();
    assert!((dt_sum - 1.0).abs() < 1e-9);
    assert_eq!(rows.last().unwrap()[4], 1.0);
    assert!(rows.iter().all(|r| r[6] <= 1.0));
    let _ = remove_dir_all(&path);
}

#[test]
fn custom_output_times() {
    // unsorted, duplicated, out of range, and missing t_end
    let path = folder("custom");
    let mut h = heat(None).unwrap();
    let mut lin: Box<dyn LinearBase> = FaerLu::new(1).unwrap();
    let mut opt = AdaptiveOptions::new(0.01, 1.0, path.to_str().unwrap()).unwrap();
    opt.out_times = vec![0.7, 0.25, -1.0, 0.0, 0.25, 2.0];
    h.phys.solve_adaptive(&mut h.vars, &mut lin, &opt).unwrap();

    let (_, rows) = read_csv(path.join("adaptive_times.csv"));
    let times: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    assert_eq!(rows.iter().map(|r| r[0]).collect::<Vec<f64>>(), vec![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(times, vec![0.0, 0.25, 0.7, 1.0]);
    let _ = remove_dir_all(&path);
}

#[test]
fn rejected_steps_keep_previous_state() {
    // a large first step is rejected; value_prev seen by the hook must always be the last accepted state
    let path = folder("reject");
    let mut h = heat(Some((0.3, 0.05))).unwrap();
    let mut lin: Box<dyn LinearBase> = FaerLu::new(1).unwrap();
    let mut opt = AdaptiveOptions::new(1.0, 1.0, path.to_str().unwrap()).unwrap();
    opt.rtol = 1e-3;
    opt.atol = 1e-6;
    opt.out_times = vec![0.3, 0.35];
    let temp = h.temp;
    let mut last = h.vars.scldom[temp].value.clone();
    let mut t_last = 0.0;
    let mut max_diff: f64 = 0.0;
    h.phys.solve_adaptive_with(&mut h.vars, &mut lin, &opt, &mut |v, t, dt| {
        let prev = &v.scldom[temp].value_prev;
        max_diff = max_diff.max(prev.iter().zip(&last).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max));
        assert!((t - dt - t_last).abs() < 1e-12);
        last = v.scldom[temp].value.clone();
        t_last = t;
    }).unwrap();
    assert_eq!(max_diff, 0.0);

    let (_, rows) = read_csv(path.join("adaptive_stats.csv"));
    assert!(rows[0][1] > 0.0);
    let dt_min = rows.iter().map(|r| r[5]).fold(f64::MAX, f64::min);
    let dt_max = rows.iter().map(|r| r[5]).fold(0.0, f64::max);
    assert!(dt_max > 10.0 * dt_min);
    let _ = remove_dir_all(&path);
}

#[test]
fn fixed_matches_solve() {
    let path = folder("fixed");
    let mut a = heat(Some((0.1, 0.05))).unwrap();
    let mut lin: Box<dyn LinearBase> = FaerLu::new(1).unwrap();
    a.phys.solve(&mut a.vars, &mut lin, 0.01, 50, 50, 30, 1e-10, 1.0).unwrap();

    let mut b = heat(Some((0.1, 0.05))).unwrap();
    let mut opt = AdaptiveOptions::new(0.01, 0.5, path.to_str().unwrap()).unwrap();
    opt.fixed = true;
    opt.tol = 1e-10;
    opt.out_times = vec![0.1, 0.15];
    b.phys.solve_adaptive(&mut b.vars, &mut lin, &opt).unwrap();

    let ta = &a.vars.scldom[a.temp].value;
    let tb = &b.vars.scldom[b.temp].value;
    let t_max = ta.iter().cloned().fold(0.0, f64::max);
    let diff = ta.iter().zip(tb).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
    assert!(t_max > 0.1);
    assert!(diff < 1e-10 * t_max);

    let (_, rows) = read_csv(path.join("adaptive_stats.csv"));
    assert_eq!(rows.len(), 50);
    let _ = remove_dir_all(&path);
}
