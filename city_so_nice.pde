class Posn {
  float x,y;
  int idx = -1;
  Posn (float _x, float _y) {
    x = _x;
    y = _y;
  }
}

class VC {
  color fill_c;
  color stroke_c;
  float stroke_w;
  
  void set_up() {
    push();
    fill(fill_c);
    stroke(stroke_c);
    strokeWeight(stroke_w);
  }
  
  void clean_up() {
    pop();
  }
}

interface Renderable {
  void display();
  VC vc=null;
}

class PQ {
  int elem[] = null;
  float wt[] = null;
  int maxsize = 0;
  int size = 0;
  
  PQ (int _maxsize) {
    elem = new int[_maxsize];
    wt = new float[_maxsize];
    maxsize = _maxsize;
  }
  
  private int _lidx (int k) {
    return 2*k + 1;
  }
  
  private int _ridx (int k) {
    return 2*k + 2;
  }
  
  private int _pidx (int k) {
    return (k - 1)/2;
  }
  
  private void _swap (int j, int k) {
      int _x = elem[j];
      float _w = wt[j];
      
      elem[j] = elem[k];
      wt[j] = wt[k];
      
      elem[k] = _x;
      wt[k] = _w;
  }
    
  private void _sift_up (int idx) {
    while (idx > 0 && wt[_pidx(idx)] > wt[idx]) {
      _swap(_pidx(idx),idx);
      idx = _pidx(idx);
    }
  }
  
  private void _sift_down (int idx) {
    int _idx = idx;
    int lidx = _lidx(idx);
    int ridx = _ridx(idx);
    
    if (lidx < size && wt[lidx] <= wt[idx]) {
      _idx = lidx;
    }
    
    if (ridx < size && wt[ridx] <= wt[idx]) {
      if (wt[ridx] <= wt[lidx]) {
        _idx = ridx;
      }
    }
    
    if (idx != _idx) {
      _swap(idx,_idx);
      _sift_down(_idx);
    }
  }
  
  void push (float _w, int _x) {
    int idx = size++;
    
    wt[idx] = _w;
    elem[idx] = _x;
    
    _sift_up(idx);
  }
  
  int pop () {
    int rslt = elem[0];
    int idx = --size;

    wt[0] = wt[idx];
    elem[0] = elem[idx];

    _sift_down(0);
    return rslt;
  }
  
  void change_w (int _elem, float w1) {
    for (int k=0; k<size; ++k) {
      /* HACK: O(n) implementation */
      if (elem[k] == _elem) {
        float w0 = wt[k];
        wt[k] = w1;
        
        if (w1 > w0) _sift_down(k);
        if (w1 < w0) _sift_up(k);
        
        break;
      }
    }
  }

  private void _write(int idx, int depth) {
    for (int k=0; k<depth; ++k) {
      print("  ");
    }
    
    print(" (", elem[idx], ": ", wt[idx]);
    if (_lidx(idx) < size) {
      print("\n");
      _write(_lidx(idx), depth+1);
    }
    if (_ridx(idx) < size) {
      print("\n");
      _write(_ridx(idx), depth+1);
    }
    print(" )");
    
    if (depth == 0) {
      print("\n");
    }
  }

  void write () {
    if (size > 0) {
      _write(0,0);
    }
  }
}

class Graph {
  private boolean adj[][];
  private boolean show[];
  Posn loc[];
  
  private int idx = 0;
  int size = 0;
  
  Graph (int n) {
    loc = new Posn[n];
    adj = new boolean[n][n];
    show = new boolean[n];
    size = n;
  }
  
  void add_to (float _x, float _y) {
    idx = ++idx % size;
    insert_at(idx,_x,_y);
  }
  
  void insert_at (int k, float _x, float _y) {
    loc[k] = new Posn(_x,_y);
    loc[k].idx = k;
    show[k] = false;
  }
  
  void unmask (int k) {
    show[k] = true;
  }
  
  void set_adj_ij (int i, int j, boolean _val) {
    adj[i][j] = _val;
    adj[j][i] = _val;
  }
  
  void display () {
    push();
    ellipseMode(CENTER);
    
    for (int i=0; i<size; ++i) {
      if (show[i]) {
        circle(loc[i].x, loc[i].y, 2.5);
        // text(str(i),loc[i].x,loc[i].y);
      }
    }
    
    for (int i=0; i<size; ++i) {
      for (int j=i+1; j<size; ++j) {
        if (adj[i][j] &&
            show[i] &&
            show[j]) {
          line(loc[i].x, loc[i].y,
               loc[j].x, loc[j].y);
        }
      }
    }
    
    pop();
  }
}

Posn get_loc (Posn ref_a, Posn ref_b, float coord) {
  Posn loc = new Posn(
    lerp(ref_a.x, ref_b.x, coord),
    lerp(ref_a.y, ref_b.y, coord));
  return loc;
}

PVector get_ortho (PVector ref, int dir) {
  PVector ortho =
    new PVector(1,-ref.x/ref.y);
  ortho.normalize();
  ortho.mult(dir);
  return ortho;
}

PVector get_displace (PVector ref_vec, float delta, int offset_dir) {
  float dx;
  float dy;

  PVector ortho = get_ortho(ref_vec, offset_dir);
  dx = delta*ortho.x;
  dy = delta*ortho.y;
  
  return new PVector(dx,dy);
}

class Location implements Renderable {
  int idx_a,idx_b;
  float offset;
  float coord;
  int dir;

  VC vc = null;
  PVector orig;
  PVector ortho;
  PVector displace;
  
  Posn loc;
  Posn ref;
  
  private float dim1=1.;
  private float dim2=1.;
  
  Location (Posn ref_a, Posn ref_b, float _coord, float _offset, int _dir) {
    idx_a = ref_a.idx;
    idx_b = ref_b.idx;
    
    offset = _offset;
    coord = _coord;
    dir = _dir;
    
    orig = new PVector(ref_b.x - ref_a.x, ref_b.y - ref_a.y);
    ortho = get_ortho(orig, dir);
    displace = get_displace(orig, offset, dir);
    
    ref = get_loc(ref_a, ref_b, coord);
    loc = new Posn(ref.x + displace.x,
                   ref.y + displace.y);
                   
    dim1 = 1. + random(7.5);
    dim2 = 1. + random(7.5);
  }
  
  void display () {
    if (vc != null) { vc.set_up(); }
    else { push(); }
    
    rectMode(CENTER);
    PVector p1 = orig;  p1.setMag(dim1);
    PVector p2 = ortho; p2.setMag(dim2);
    
    PVector q1 = new PVector(loc.x + p1.x + p2.x, loc.y + p1.y + p2.y),
            q2 = new PVector(loc.x + p1.x - p2.x, loc.y + p1.y - p2.y),
            q4 = new PVector(loc.x - p1.x + p2.x, loc.y - p1.y + p2.y),
            q3 = new PVector(loc.x - p1.x - p2.x, loc.y - p1.y - p2.y);

    quad(q1.x, q1.y,
         q2.x, q2.y,
         q3.x, q3.y,
         q4.x, q4.y);
    
    if (vc != null) { vc.clean_up(); }
    else { pop(); }
  }
}


int lidx_to_nidx (Map m, int lidx, Posn _ref) {
  Location loc = m.loc.get(lidx);
  PVector ref = new PVector(_ref.x, _ref.y);
  /* doesn't work for the random map!!! */
  
  Posn _a = m.g.loc[loc.idx_a];
  Posn _b = m.g.loc[loc.idx_b];
  
  PVector a = new PVector(_a.x,_a.y);
  PVector b = new PVector(_b.x,_b.y);
  
  float da = a.dist(ref);
  float db = b.dist(ref);
  
  return db >= da ? loc.idx_a : loc.idx_b;
}

boolean same_segment(Location la, Location lb) {
  boolean min_eq = min(la.idx_a,la.idx_b) == min(lb.idx_a,lb.idx_b);
  boolean max_eq = max(la.idx_a,la.idx_b) == max(lb.idx_a,lb.idx_b);
  return min_eq && max_eq;
}


ArrayList<Posn> find_path (Map m, int srcl, int dstl) {
  ArrayList<Posn> path = new ArrayList<Posn>();
  
  Posn srcp = m.loc.get(srcl).ref;
  Posn dstp = m.loc.get(dstl).ref;
  
  path.add(srcp);
  
  int srcn = lidx_to_nidx(m, srcl, dstp);
  int dstn = lidx_to_nidx(m, dstl, srcp);
  
  if (! same_segment(m.loc.get(srcl),m.loc.get(dstl))) {
    IntList np = djikstra(m.g, srcn, dstn);
    for (int nk : np) {
      path.add(m.g.loc[nk]);
    }
  }
  
  path.add(dstp);
  return path;
}


IntList djikstra (Graph g, int src, int dst) {
  IntList p = new IntList();
  IntList q = new IntList();
  
  int par[] = new int[g.size];
  float memo[] = new float[g.size];
  
  for (int k=0; k<g.size; ++k) {
    memo[k] = -1.0;
    par[k] = -1;
  }

  PQ minq = new PQ(g.size);
  minq.push(0,src);
  memo[src] = 0;

  int curr = src;
  while (curr != dst && minq.size > 0) {
    PVector x0 = new PVector(
      g.loc[curr].x, g.loc[curr].y);
    
    for (int k=0; k<g.size; ++k) {
      if (g.adj[curr][k]) {
        PVector x1 = new PVector(
          g.loc[k].x, g.loc[k].y);
        
        float dist =
          x0.dist(x1) +
          memo[curr];
        
        if (memo[k] < 0) {
          minq.push(dist,k);
          memo[k] = dist;
          par[k] = curr;
        } else if (dist <= memo[k]) {
          minq.change_w(k,dist);
          memo[k] = dist;
          par[k] = curr;
        }
      }
    }
    
    curr = minq.pop();
  }
  
  if (curr == dst) {
    while (curr != src) {
      q.append(curr);
      curr = par[curr];
    }
    q.append(src);
  }
  
  for (int k=q.size()-1; k>=0; --k) {
    p.append(q.get(k));
  }
  return p;
}


class Map {
  Graph g;
  ArrayList<Location> loc;
  
  Map (Graph _g) {
    loc = new ArrayList<Location>();
    g = _g;
  }

  void display () {
    g.display();
    //int k=0;
    for (Location _loc : loc) {
      _loc.display();
      //text(str(k++),_loc.loc.x,_loc.loc.y-5);
    }
  }
}

class Automata implements Renderable {
  ArrayList<Posn> q = new ArrayList<Posn>();
  boolean active=false;
  Posn src, dst;
  int state = -1;
  
  float coord=0;
  float dxdt=0;
  VC vc=null;
  
  private PVector dir=null;
  private PVector loc=null;
  
  /* HACK HACK HACK */
  private float displacement_factor=0.0;
  
  Automata (Posn p0, float dxdt0) {
    src = p0;
    dst = p0;
    dir = new PVector (0,0);
    coord = 1.0;
    dxdt = dxdt0;
    loc = new PVector(p0.x,p0.y);
  }
  
  private void _next() {
    if (q.size() > 0) {
      coord = 0.0;
        
      src = dst;
      dst = q.get(0);
      q.remove(0);
        
      dir = new PVector(dst.x-src.x,
                        dst.y-src.y);
      dir.normalize();
    } else {
      coord = 1.0;
      dir = new PVector(0,0);
      active = false;
    }
  }
  
  void step_by (float px) {
    if (active) {
      PVector delta = dir;
      delta.setMag(px);
      loc = loc.add(delta);
      coord = new PVector(loc.x-dst.x,
                          loc.y-dst.y).dot(dir);
      if (coord >= 0) {
        _next();
      }
    }
  }
  
  void step (float dt) {
    if (active) {
      coord += dxdt*dt;
      if (coord > 1.0) {
        _next();
      }
      loc = new PVector(
        lerp(src.x,dst.x,coord),
        lerp(src.y,dst.y,coord));
    }
  }
  
  void send_to (Posn p) {
    active = true;
    q.add(p);
  }
  
  void set_at (Posn p, int st) {
    active = false;
    q.clear();

    coord = 1.0;
    src = p;
    dst = p;
    state = st;

    dir = new PVector(0,0);
    loc = new PVector(p.x,p.y);
  }
  
  void display () {
    if (vc != null) { vc.set_up(); }
    else { push(); }

    /* HACK HACK HACK */
    PVector displ = new PVector(0,0);
    float _adf = abs(displacement_factor); 
    if (_adf > 0) {
      int _sdf = int(displacement_factor/_adf);
      PVector _delta = get_displace(dir,_adf,_sdf);
      displ.x += _delta.x;
      displ.y += _delta.y;
    }

    ellipseMode(CENTER);
    circle(loc.x+displ.x, loc.y+displ.y, 1.5);
    
    // text(str(state), loc.x+10, loc.y);
    if (vc != null) { vc.clean_up(); }
    else { pop(); }
  }
}




/*
 * run calls and utility functions
 */

float _rcoord(int dim, float pad) {
  float retval;
  int p = int(dim*pad);
  retval = p + random(dim - 2*p);
  return retval;
}

class TestMap extends Map {
  /* intended for DIM x DIM test canvas */
  TestMap (Graph _g, int n, int dim) {
    super(_g);
    
    for (int i=0; i<n; ++i) {
      if (0 != i) g.set_adj_ij(i-1,i,true);
      g.insert_at(i,_rcoord(dim,0.1),_rcoord(dim,0.1));
      g.unmask(i);
    }
    g.set_adj_ij(n-1,0,true);
    
    for (int i=0; i<n; ++i) {
      int j = (i == n-1) ? 0 : (i+1);
      
      Location l1 = new Location(g.loc[i],g.loc[j],random(0.1,0.9),5,-1);
      Location l2 = new Location(g.loc[i],g.loc[j],random(0.1,0.9),5, 1);
      
      loc.add(l1);
      loc.add(l2);
    }
  }
}



class GridMap extends Map {
  int nrow=0; float rsep=0;
  int ncol=0; float csep=0;
  
  float pad=0;
  
  GridMap (Graph _g, int _nrow, int _ncol, int nloc,
           int row_sep, int col_sep,
           float noise_parm,
           float pad_by) {
    super(_g);
    
    nrow=_nrow; rsep=row_sep;
    ncol=_ncol; csep=col_sep;
    
    float _w = pad_by * (ncol-1)*col_sep;
    float _h = pad_by * (nrow-1)*row_sep;
    
    float xloc=_w,yloc=_h;
    for (int j=0; j<nrow; ++j) {
      yloc += row_sep;
      xloc=_w;
      
      for (int k=0; k<ncol; ++k) {
         xloc += col_sep;
         
         int idx = _rc_to_idx(j,k);
         g.loc[idx] = new Posn(xloc,yloc);
         g.loc[idx].idx = idx;
         g.unmask(idx);
         
         if (k > 0) {
           g.set_adj_ij(idx, _rc_to_idx(j,k-1), true);
         }
         
         if (k < ncol-1) {
           g.set_adj_ij(idx, _rc_to_idx(j,k+1), true);
         }
         
         if (j > 0) {
           g.set_adj_ij(idx, _rc_to_idx(j-1,k), true);
         }
         
         if (j < ncol-1) {
           g.set_adj_ij(idx, _rc_to_idx(j+1,k), true);
         }
      }
    }
    
    for (int j=0; j<nrow; ++j) {
      for (int k=0; k<ncol; ++k) {
        int idx = _rc_to_idx(j,k);
        g.loc[idx].x += 2*(noise(random(100), random(100))-0.5)*noise_parm*row_sep;
        g.loc[idx].y += 2*(noise(random(100), random(100))-0.5)*noise_parm*col_sep;
      }
    }
    
    for (int j=0; j<nloc; ++j) {
      int ref_a = int(random(nrow*ncol));
      int ref_b = -1;

      while (ref_b < 0) {
        float p = random(1);
        int ridx = int(ref_a/ncol);
        int cidx = ref_a % ncol;
        
        if (p <= 0.25) {
          if (cidx != 0) {
            ref_b = _rc_to_idx(ridx,cidx-1);
          }
        } else if (p <= 0.5) {
          if (cidx != ncol-1) {
            ref_b = _rc_to_idx(ridx,cidx+1);
          }
        } else if (p <= 0.75) {
          if (ridx != 0) {
            ref_b = _rc_to_idx(ridx-1,cidx);
          }
        } else {
          if (ridx != nrow-1) {
            ref_b = _rc_to_idx(ridx+1,cidx);
          }
        }
      }

      int dir = (random(1) <= 0.5) ? -1 : 1;
      Location l = new Location(g.loc[ref_a],
                                g.loc[ref_b],
                                random(0.1,0.9),10,
                                dir);      
      loc.add(l);
    }
  }
  
  private int _rc_to_idx(int r, int c) {
    return r*ncol + c;
  }
}

void test_pq() {
  // 4 1 8 7 6 5 9 2 3 
  PQ pq = new PQ(100);
  pq.push(10.5,  9);
  pq.push( 8.25, 8);
  pq.push(11.75, 7);
  pq.push( 9.95, 6);
  pq.push( 8.35, 5);
  pq.push( 0.5,  4);
  pq.push(90.75, 3);
  pq.push(11.25, 2);
  pq.push(6.5,   1);
  pq.change_w(5, 10.05);
  pq.change_w(7, 9.25);
  while (pq.size > 0) {
    println(pq.pop());
  }
}

Map m;
Automata a[];
int nauto = 1500; // 250;
int nloc = 350; // 50;
char lbl[] = new char[nloc];
int clock = 0;
int start_tone=225;

color rand_c(int minval) {
  int ub = 255-minval;
  int r = int(random(ub)),
      g = int(random(ub)),
      b = int(random(ub));
  return color(r+minval,g+minval,b+minval);
}

void setup() {
  surface.setTitle("The city so nice they named it twice");
  // size(500,500);
  size(650,600);
  // int dim = 500;
  a = new Automata[nauto];
  
  // m = new TestMap(new Graph(nelem), nelem, dim);
  m = new GridMap(new Graph(70),
                  14, 
                  5,
                  nloc,
                  35, 100,
                  0.25, 0.05);
  
  /*
  for (int k=0; k<nauto; ++k) {
    int st = int(random(nelem+1));
    a[k] = new Automata(m.g.loc[st],0);
    a[k].set_at(m.g.loc[st],st);
  }
  */
  
  for (int k=0; k<nauto; ++k) {
    int st = int(random(m.loc.size()));
    a[k] = new Automata(m.loc.get(st).ref,0);
    a[k].set_at(m.loc.get(st).ref,st);
  }
  
  for (Automata _a : a) {
    // color col = rand_c(85);
    color col = color(start_tone);
    _a.vc = new VC();
    _a.vc.fill_c=col;
    _a.vc.stroke_c=col;
    _a.vc.stroke_w=2.0;
  }
  
  for (Location loc : m.loc) {
    loc.vc = new VC();
    //color col = color(start_tone);
    color col = rand_c(85);
    loc.vc.fill_c = col;
    loc.vc.stroke_c = col;
    loc.vc.stroke_w = 2.0;
  }
  
  String lblstr = new
    String("AzByCxDwEvFuGtHsIrJqKpLoMnNmOlPlQjRiShTgUfVeWdXcYbZa");
  for (int k=0; k<nloc; ++k) {
    lbl[k] = lblstr.charAt(k%52);
  }
}

void draw() {
  background(45);
  stroke(125);

  clock = (clock+1)%1000;
  float prob_t = 0.0125*(cos(clock/250. * TWO_PI)+1)/2;

  /*
  for (int k=0; k<nauto; ++k) {
    if (! a[k].active) {
      if (random(1) <= 0.0025) {
        int st = int(random(m.g.size));
        if (st != a[k].state) {
          IntList p = djikstra(m.g,a[k].state,st);
          for (int x : p) {
            a[k].send_to(m.g.loc[x]);
          }
          a[k].state = st;
        }
      }
    }
  }
  */

  /*
  int ctr[] = new int[nloc];
  color tmp[] = new color[nloc];
  for (int k=0; k<nloc; ++k) {
    tmp[k]=color(start_tone);
    ctr[k]=0;
  }
  */
  
  for (int k=0; k<nauto; ++k) {
    if (! a[k].active) {
      a[k].displacement_factor = 0;
      if (random(1) <= prob_t) {
        int st = int(random(m.loc.size()));
        if (st != a[k].state) {
          ArrayList<Posn> pt = find_path(m,a[k].state,st);
          for (Posn p : pt) {
            a[k].send_to(p);
          }
          a[k].state = st;
          a[k].dxdt = (random(2)+1)/2.5;
          
          a[k].vc.stroke_c = m.loc.get(st).vc.stroke_c;
          a[k].vc.fill_c = m.loc.get(st).vc.fill_c;
        }
      } else {
        /*
          int st = a[k].state;
          assert(a[k].vc.stroke_c ==
                 a[k].vc.fill_c);
          tmp[st] = lerpColor(tmp[st],a[k].vc.fill_c,1/float(++ctr[st]));
         */
      }
    }
  }
  
  for (int k=0; k<nauto; ++k) {
    if (a[k].active) {
      a[k].vc.set_up();
      /*
      text(str(lbl[a[k].state]),
               a[k].loc.x+10,
               a[k].loc.y);
      */
      a[k].step_by(a[k].dxdt);
      a[k].display();
      a[k].vc.clean_up();
    }
  }
  
  for (int k=0; k<m.loc.size(); ++k) {
    /*
    m.loc.get(k).vc.stroke_c = tmp[k];
    m.loc.get(k).vc.fill_c = tmp[k];
    */

    /*
    m.loc.get(k).vc.set_up();
    text(str(ctr[k]),
         m.loc.get(k).loc.x+5,
         m.loc.get(k).loc.y+5);
    m.loc.get(k).vc.clean_up();
    */
  }
  
  m.display();
}
