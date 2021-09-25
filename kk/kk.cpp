#include <Core/Core.h>

using namespace Upp;


class Parent {
public:
	void fun() {data--;}
	
	void kk();
	
protected:
	int64 data = 0;

friend class SpecializationA;
friend class SpecializationB;
//friend class SpecializationC;		No necesario
};

class Parent_Format {
public:
	virtual void myfunction() = 0;
};

class SpecializationA : public Parent_Format {   
public:
	SpecializationA(Parent &inter) : inter(inter) {}
    virtual void myfunction() {inter.data--;}

private:
	Parent &inter;
};

class SpecializationB {
public:
    static void myfunction(Parent &d) {d.data--;}; 
};

class SpecializationC : public Parent {
public:
	void myfunction() {data--;}
};

void Parent::kk() {
	static_cast<SpecializationC&>(*this).myfunction(); 
}
	
CONSOLE_APP_MAIN
{
	Parent a;
	
	a.fun();
	TimeStop t;
	int num = 90000000;
	
	t.Reset();
	for (int i = 0; i < num; ++i)
		a.fun();
	Cout() << "\n" << t.Seconds();
	
	t.Reset();
	for (int i = 0; i < num; ++i)
		SpecializationA(a).myfunction();
	Cout() << "\n" << t.Seconds();
	
	t.Reset();
	for (int i = 0; i < num; ++i)
		SpecializationB::myfunction(a);
	Cout() << "\n" << t.Seconds();
	
	t.Reset();
	for (int i = 0; i < num; ++i)
		static_cast<SpecializationC&>(a).myfunction();
	Cout() << "\n" << t.Seconds();
	
	ReadStdIn();
}
