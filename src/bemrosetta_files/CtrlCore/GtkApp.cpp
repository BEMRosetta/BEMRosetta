#include <CtrlCore/CtrlCore.h>

#ifdef GUI_GTK

#define CATCH_ERRORS 0

namespace Upp {

#define LLOG(x) // DLOG(x)

#if CATCH_ERRORS
void CatchError(const gchar *log_domain,
             GLogLevelFlags log_level,
             const gchar *message,
             gpointer user_data)
{
	RLOG((const char *)message);
	__BREAK__;
}
#endif

void _DBG_Ungrab(void)
{   // This is a special nasty hack to make possible to ungrab mouse by debugger (see ide/Debuggers/PrettyPrinters.py)
	gdk_pointer_ungrab(GDK_CURRENT_TIME);
}

void Ctrl::PanicMsgBox(const char *title, const char *text)
{
	LLOG("PanicMsgBox " << title << ": " << text);
	if(gdk_pointer_is_grabbed())
		gdk_pointer_ungrab(CurrentTime);
	char m[2000];
	*m = 0;
	if(system("which gxmessage") == 0)
		strcpy(m, "gxmessage -center \"");
	else
	if(system("which kdialog") == 0)
		strcpy(m, "kdialog --error \"");
	else
	if(system("which xmessage") == 0)
		strcpy(m, "xmessage -center \"");

	if(*m) {
		strcat(m, title);
		strcat(m, "\n");
		strcat(m, text);
		strcat(m, "\"");
		IGNORE_RESULT(system(m));
	}
	else {
		_DBG_Ungrab();
		GtkWidget *dialog = gtk_message_dialog_new(NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR,
		                                           GTK_BUTTONS_CLOSE, "%s: %s", title, text);
		gtk_dialog_run(GTK_DIALOG (dialog));
		gtk_widget_destroy(dialog);
	}
	__BREAK__;
}

void InitGtkApp(int argc, char **argv, const char **envptr)
{
	LLOG(rmsecs() << " InitGtkApp");
#ifdef _MULTITHREADED
#if !GLIB_CHECK_VERSION(2, 32, 0)
    if(!g_thread_supported())
        g_thread_init(NULL);
#endif
	gdk_threads_set_lock_functions(EnterGuiMutex, LeaveGuiMutex);
	gdk_threads_init();
	EnterGuiMutex();
#endif
	gtk_init(&argc, &argv);
	Ctrl::GlobalBackBuffer();
	Ctrl::ReSkin();
	g_timeout_add(20, (GSourceFunc) Ctrl::TimeHandler, NULL);
	InstallPanicMessageBox(Ctrl::PanicMsgBox);
	gdk_window_add_filter(NULL, Ctrl::RootKeyFilter, NULL);
#if CATCH_ERRORS
	g_log_set_default_handler (CatchError, 0);
#endif
}

void ExitGtkApp()
{
#ifdef _MULTITHREADED
	LeaveGuiMutex();
#endif
}

}

#endif
