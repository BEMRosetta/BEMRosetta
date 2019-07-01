//$ class Upp::TopWindow : Ctrl {
	ImageGdk gdk_icon;
	ImageGdk gdk_largeicon;
	bool     topmost;
	
	enum { FULLSCREEN = 99 };

	void     CenterRect(Ctrl *owner);
	void     SetMode(int mode);
	void     SyncTopMost();

	static gboolean StateEvent(GtkWidget *widget, GdkEventWindowState *event, gpointer user_data);

	friend class Ctrl;
//$ };
