#ifndef _74c75f6d_1c50_4b5c_96aa_bfe42eaf6e36
#define _74c75f6d_1c50_4b5c_96aa_bfe42eaf6e36

#include "bittree_prelude.hxx"
#include "bittree_mem_defs.hxx"

namespace BitTree {
  template<class T> class Ref;

  class Ref_ {
  protected:
    struct Obj {
      size_t refs;
      void *ptr;
      void (*dtor)(Obj*);
      static void blank_dtor(Obj *me) {
        me->~Obj();
        free(me);
      }
      template<class T>
      static void inst_dtor(Obj *me) {
        reinterpret_cast<T*>(me->ptr)->~T();
        me->~Obj();
        free(me);
      }
      Obj(void *ptr, void(*dtor)(Obj*)):
        refs(1),
        ptr(ptr),
        dtor(dtor) {
      }
    };
    struct Link: Obj {
      Obj *inst;
      static void dtor(Obj *o) {
        Link *me = static_cast<Link*>(o);
        if(--me->inst->refs == 0)
          me->inst->dtor(me->inst);
        me->~Link();
        free(me);
      }
      Link(void *ptr, Obj *inst): Obj(ptr, dtor), inst(inst) {}
    };
    struct Array: Obj {
      size_t n;
      template<class T>
      static void dtor(Obj *o) {
        Array *me = static_cast<Array*>(o);
        for(size_t i=0; i < me->n; i++)
          (reinterpret_cast<T*>(me->ptr) + i)->~T();
        me->~Array();
        free(me);
      }
      Array(void *ptr, void(*dtor)(Obj*), size_t n): Obj(ptr, dtor), n(n) {}
    };
    
    template<class pA, class pB>
    struct implicit_caster {
      static void* cast(void *a) { pB b = (pA)a; return (void*)b; }
    };
    template<class pA, class pB>
    struct const_caster {
      static void* cast(void *a) { return (void*)const_cast<pB>((pA)a); }
    };
    template<class pA, class pB>
    struct static_caster {
      static void* cast(void *a) { return (void*)static_cast<pB>((pA)a); }
    };
    template<class pA, class pB>
    struct dynamic_caster {
      static void* cast(void *a) { return (void*)dynamic_cast<pB>((pA)a); }
    };
    template<class pA, class pB>
    struct reinterpret_caster {
      static void* cast(void *a) { return (void*)reinterpret_cast<pB>((pA)a); }
    };
    template<class caster>
    static Obj* acquire(Obj *a) {
      if(a == 0x0) return 0x0;
      void *b_ptr = caster::cast(a->ptr);
      if(a->ptr == b_ptr) {
        a->refs += 1;
        return a;
      }
      Obj *inst = a->dtor == Link::dtor
        ? static_cast<Link*>(a)->inst
        : a;
      inst->refs += 1;
      return new Link(b_ptr, inst);
    }
    template<class T>
    static Ref<T> make(Obj *o) {
      return Ref<T>(o);
    }
  protected:	
    Obj *o;
    Ref_(Obj *o): o(o) {}
  protected:
    template<class pA, class pB>
    void assign(Obj *o1) {
      if(o != o1) {
        if(o && 0 == --o->refs)
          o->dtor(o);
        o = acquire<implicit_caster<pA,pB> >(o1);
      }
    }
  public:
    Ref_(): o(0x0) {}
    Ref_(const Ref_ &that):
      o(acquire<implicit_caster<void*,void*> >(that.o)) {
    }
    ~Ref_() {
      if(o && 0 == --o->refs) 
        o->dtor(o);
    }
    Ref_& operator=(const Ref_ &that) {
      this->assign<void*,void*>(that.o);
      return *this;
    }
    void nullify() {
      if(o) {
        if(0 == --o->refs)
          o->dtor(o);
        o = 0x0;
      }
    }
    template<class U>
    Ref<U> cast_reinterpret() const {
      return Ref<U>(acquire<reinterpret_caster<void*,U*> >(o));
    }
    void* alloc(const Mem &m) {
      Mem m1 = Mem::of<Obj>();
      size_t off = m1.concat(m);
      void *inst = malloc(m1.size);
      void *ptr = reinterpret_cast<char*>(inst) + off;
      this->nullify();
      this->o = new(inst) Obj(ptr, Obj::blank_dtor);
      return ptr;
    }
    template<class T>
    Ref<T> promote() {
      Obj *o = this->o;
      if(o) {
        this->o = 0x0;
        o->dtor = Obj::inst_dtor<T>;
      }
      return Ref<T>(o);
    }
    
    template<class T>
    static Ref<T> new1() {
      Ref<T> ref;
      new(ref.alloc()) T();
      return ref;
    }
    template<class T, class X0>
    static Ref<T> new1(X0 x0) {
      Ref<T> ref;
      new(ref.alloc()) T(x0);
      return ref;
    }
    template<class T>
    static Ref<T> new_array(size_t n) {
      Mem m = Mem::of<Array>();
      size_t off = m.concat(Mem::array(Mem::of<T>(), n));
      void *obj = malloc(m.size);
      void *ptr = reinterpret_cast<char*>(obj) + off;
      for(size_t i=0; i < n; i++)
        new(reinterpret_cast<char*>(ptr) + i*sizeof(T)) T();
      return Ref<T>(new(obj) Array(ptr, Array::dtor<T>, n));
    }
  };

  template<class T>
  class Ref: public Ref_ {
    friend class Ref_;
    template<class U> friend class Ref;
    Ref(Obj *o): Ref_(o) {}
  public:
    Ref(): Ref_() { }
    template<class U>
    Ref(const Ref<U> &that): Ref_(acquire<implicit_caster<U*,T*> >(that.o)) {}
    template<class U>
    Ref<T>& operator=(const Ref<U> &that) {
      this->assign<U*,T*>(that.o);
      return *this;
    }
    
    operator T*() const { return reinterpret_cast<T*>(o ? o->ptr : 0x0); }
    T* operator->() const { return reinterpret_cast<T*>(o->ptr); }
    T& operator*() const { return *reinterpret_cast<T*>(o->ptr); }
    T& operator[](size_t i) const { return reinterpret_cast<T*>(o->ptr)[i]; }
    
    template<class U>
    Ref<U> cast_const() const {
      return Ref<U>(acquire<const_caster<T*,U*> >(o));
    }
    template<class U>
    Ref<U> cast_static() const {
      return Ref<U>(acquire<static_caster<T*,U*> >(o));
    }
    template<class U>
    Ref<U> cast_dynamic() const {
      return Ref<U>(acquire<dynamic_caster<T*,U*> >(o));
    }
  private:
    Ref_::promote;
  public:
    void* alloc(Mem m = Mem::of<T>()) {
      void *x = Ref_::alloc(m);
      this->o->dtor = Obj::inst_dtor<T>;
      return x;
    }
    
    void* to_naked() {
      return o;
    }
    static Ref<T> from_naked(void *o) {
      Obj *obj = reinterpret_cast<Obj*>(o);
      obj->refs += 1;
      return Ref<T>(obj);
    }
  };
}
#endif
