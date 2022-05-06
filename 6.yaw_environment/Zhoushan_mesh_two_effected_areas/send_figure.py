import time

from PIL import ImageGrab
import smtplib
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.header import Header

endDate=time.strftime('%Y-%m-%d',time.localtime(time.time()))
# img = ImageGrab.grab()
# img.save('E:\\12.png')
my_sender = '123@qq.com' # 发件人邮箱账号
my_pass = 'tmugmde3333ad' # 发件人邮箱密码
my_user = '123@qq.com' # 收件人邮箱账号，我这边发送给自己
sender = “123@qq.com”
receivers = ['123@qq.com'] # 接收邮件，可设置为你的QQ邮箱或者其他邮箱
msgRoot = MIMEMultipart('related')
msgRoot['From'] = Header(str(endDate)+" 大盘趋势", 'utf-8')
msgRoot['To'] = Header("测试", 'utf-8')
subject = str(endDate)+' 趋势'
msgRoot['Subject'] = Header(subject, 'utf-8')

msgAlternative = MIMEMultipart('alternative')
msgRoot.attach(msgAlternative)

mail_msg = """
<p>大盘趋势</p>
<p>图片演示：</p>
<p><img src="cid:image1"></p>
"""
msgAlternative.attach(MIMEText(mail_msg, 'html', 'utf-8'))

# 指定图片为当前目录
fp = open('E:\\12.png', 'rb')
msgImage = MIMEImage(fp.read())
fp.close()

# 定义图片 ID，在 HTML 文本中引用
msgImage.add_header('Content-ID', '<image1>')
msgRoot.attach(msgImage)

try:
    smtpObj = smtplib.SMTP()
    smtpObj.connect('smtp.qq.com', 25) # 25 为 SMTP 端口号
    smtpObj.login(my_user,my_pass)
    smtpObj.sendmail(sender, receivers, msgRoot.as_string())
    print("邮件发送成功")
except smtplib.SMTPException:
    print("Error: 无法发送邮件")
